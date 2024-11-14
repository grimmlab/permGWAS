import numpy as np
import torch
import time
import scipy.stats as stats

from . import _base_model
from preprocess import data_loader
from optimize import brent


class LMM(_base_model.BaseModel):

    def __init__(self, dataset: data_loader.Dataset, batch_size: int, device: torch.device, perm: int = None,
                 perm_batch_size: int = None):
        super().__init__(dataset=dataset, batch_size=batch_size, device=device, perm=perm,
                         perm_batch_size=perm_batch_size)
        self.D = None  # eigenvalues of K
        self.U = None  # unitary matrix of eigenvectors of K
        self.freedom_deg = None  # adjusted degrees of freedom = n_samples - degrees of freedom = int
        self.Uy = None  # y after linear transformation with eigenvectors
        self.UZ = None  # fixed effects after linear transformation with eigenvectors

    def gwas(self):
        """
        Perform batch-wise computation of univariate test with linear mixed model (EMMAX):
        (1) compute spectral decomposition of K=UDU'
        (2) transform data: U'y, U'Z
        (3) estimate delta and compute variance components
        (4) calculate residual sum of squares of null model
        (5) batch-wise:
            (a) linearly transform marker
            (b) calculate effect size, residual sum of squares and standard error
            (c) calculate test statistic
        (6) calculate p-values
        Dataset:
            X: genotype matrix of shape (n,m) or (n,b) if batch-wise
            y: phenotype vector of shape (n)
            K: kinship matrix of shape (n,n)
            fixed: vector/matrix of fixed effects of shape (n) or (n,c)
        """
        start = time.time()
        self.freedom_deg = self.dataset.n_samples - self.dataset.fixed.shape[1]
        # get spectral decomposition
        self.D, self.U = torch.linalg.eigh(self.dataset.K)
        # linearly transform data, i.e. compute U'y and U'Z for fixed effects Z
        self.Uy = self.transform_input(X=self.dataset.y, U=self.U)  # shape (n)
        self.UZ = self.transform_input(X=self.dataset.fixed, U=self.U)  # shape (n) or (n,c)
        # estimate delta and compute variance components
        self.delta = self.estimate_delta(gridlen=100, logdelta_min=-10, logdelta_max=10, reml=True)
        D = self.D + self.delta
        ZD = self._zd(UZ=self.UZ, D=D)
        ZDZ = self._zdz(UZ=self.UZ, ZD=ZD)
        self.v_g, self.v_e = self.compute_var_components(D=D, UZ=self.UZ, ZD=ZD, ZDZ=ZDZ, reml=True)
        # calculate residual sum of squares of null model
        RSS_0 = self.get_rss_h0()  # shape: (1)
        self.freedom_deg -= 1
        # in batches:
        SE = []
        effect_size = []
        test_stat = []
        for batch in range(int(np.ceil(self.dataset.n_snps / self.batch_size))):
            # set bounds for SNP batch
            lower_bound, upper_bound = self._bounds(batch_size=self.batch_size, batch=batch)
            # load and transform batch of SNPs
            US = self._s_matrix(lower_bound=lower_bound, upper_bound=upper_bound)  # shape: (n,b)
            # transform data
            US = self.transform_input(X=US, U=self.U)
            # calculate effect size, residual sum of squares and standard error
            RSS_1, stds, betas = self.get_rss_and_se(D=D, S=US, ZD=ZD, ZDZ=ZDZ)
            SE.append(stds.to(torch.device("cpu")))
            effect_size.append(betas.to(torch.device("cpu")))
            # calculate test statistic
            test_stat.append(self.get_f_score(rss0=RSS_0, rss1=RSS_1).to(torch.device("cpu")))
            # free GPU space
            if self.device.type != "cpu":
                with torch.cuda.device(self.device):
                    del RSS_1
                    del US
                    del stds
                    del betas
                    torch.cuda.empty_cache()
        self.SE = torch.cat(SE, dim=0)  # shape: (m)
        self.effect_size = torch.cat(effect_size, dim=0)  # shape: (m)
        self.test_stat = torch.cat(test_stat, dim=0)  # shape: (m)
        time_test_stats = time.time()
        print("Have test statistics of %d SNPs. Elapsed time: %f" % (self.test_stat.shape[0], time_test_stats - start))
        print("Calculate P-values now")
        # compute p-values
        self.p_value = torch.tensor(list(map(self.get_p_value, self.test_stat)))
        print("Have P-values. Elapsed time: ", time.time() - time_test_stats)
        if self.device.type != "cpu":
            with torch.cuda.device(self.device):
                del D
                del ZD
                del ZDZ
                del self.dataset.K
                torch.cuda.empty_cache()

    def perm_gwas(self, perm_method: str = 'x', adj_p_value: bool = False):
        """
        Perform batch-wise computation of permutation-based test with linear mixed model (EMMAX):
        reuse spectral decomposition of K=UDU'
        perm method y:
            (1) permute phenotype p times
            (2) transform data: U'y
            (3) estimate delta and compute variance components for each permutation
            (4) calculate residual sum of squares of null model
            (5) batch-wise:
                (a) linearly transform marker
                (b) calculate residual sum of squares
                (c) calculate test statistic
        perm method x:
            (1) permute fixed effects p times
            (2) transform data: U'Z
            (3) estimate delta and compute variance components for each permutation
            (4) calculate residual sum of squares of null model
            (5) batch-wise:
                (a) permute marker p times
                (b) linearly transform marker
                (c) calculate residual sum of squares
                (d) calculate test statistic
        (6) calculate minimal p-values for Westfall-Young permutation-based threshold
        optional: (7) calculate permutation-based p-values
        Dataset:
            X: genotype matrix of shape (n,m) or (n,b) if batch-wise
            y: phenotype vector of shape (n)
            K: kinship matrix of shape (n,n)
            fixed: vector/matrix of fixed effects of shape (n) or (n,c)

        :param perm_method: y to permute phenotype or x to permute fixed effects + marker
        :param adj_p_value: if True compute adjusted p-values, default is false
        """
        start = time.time()
        if self.test_stat is None:
            raise Exception('Need to first calculate true test statistics using LMM.gwas().')
        self.freedom_deg = self.dataset.n_samples - self.dataset.fixed.shape[1]
        self.seeds = self.perm_seeds()
        if perm_method == 'y':
            # compute permutations of y
            self.Uy = self.permute(data=self.dataset.y)  # shape: (n,p)
            self.Uy = torch.unsqueeze(torch.t(self.transform_input(X=self.Uy, U=self.U)), 2)  # shape: (p,n,1)
            # estimate variance components for each permutation
            self.delta = self.estimate_delta_perm(gridlen=100, logdelta_min=-10, logdelta_max=10, reml=True)
            self.D = self._d_delta(delta=self.delta, batch_size=self.perm)  # shape: (p,1,n)
            self.UZ = self.get_3d_copy(v=self.UZ, batch_size=self.perm)  # shape: (p,n,c)
            ZD = self._zd(UZ=self.UZ, D=self.D)  # shape: (p,c,n)
            ZDZ = self._zdz(UZ=self.UZ, ZD=ZD)  # shape: (p,c,c)
            v_g, _ = self.compute_var_components(D=self.D, UZ=self.UZ, ZD=ZD, ZDZ=ZDZ, reml=True)  # shape: (p)
        elif perm_method == 'x':
            self.Uy = self.get_3d_copy(v=self.Uy, batch_size=self.perm)  # shape: (p,n,1)
            if self.dataset.fixed.shape[1] > 1:
                # permute and transform fixed effects
                self.UZ = self.permute(data=self.dataset.fixed)  # shape: (p,n,c)
                self.UZ = self.transform_input(X=self.UZ, U=self.U)  # shape: (p,n,c)
                # estimate variance components
                self.delta = self.estimate_delta_perm(gridlen=100, logdelta_min=-10, logdelta_max=10, reml=True)
                self.D = self._d_delta(delta=self.delta, batch_size=self.perm)  # shape: (p,1,n)
                ZD = self._zd(UZ=self.UZ, D=self.D)  # shape: (p,c,n)
                ZDZ = self._zdz(UZ=self.UZ, ZD=ZD)  # shape: (p,c,c)
                v_g, _ = self.compute_var_components(D=self.D, UZ=self.UZ, ZD=ZD, ZDZ=ZDZ, reml=True)  # shape: (p)
            else:
                # reuse UZ, delta, sigma and get 3D copies
                self.D = self._d_delta(delta=self.delta, batch_size=self.perm)  # shape: (p,1,n)
                self.UZ = self.get_3d_copy(v=self.UZ, batch_size=self.perm)  # shape: (p,n,c)
                ZD = self._zd(UZ=self.UZ, D=self.D)  # shape: (p,c,n)
                ZDZ = self._zdz(UZ=self.UZ, ZD=ZD)  # shape: (p,c,c)
                v_g = self.v_g.repeat(self.perm)  # shape: (p)
        else:
            raise Exception('Choose either permutation method x or y.')
        # calculate rss for null model
        RSS_0 = self.get_rss_h0().repeat(self.perm)  # shape: (p)
        self.freedom_deg -= 1
        if self.device.type != "cpu":
            with torch.cuda.device(self.device):
                del self.delta
                del self.dataset.y
                del self.dataset.fixed
                torch.cuda.empty_cache()
        var_comp_time = time.time()
        print("Have variance components. Elapsed time: ", var_comp_time - start)
        test_stat = []
        for batch in range(int(np.ceil(self.dataset.n_snps / self.perm_batch_size))):
            # set bounds for SNP batch
            lower_bound, upper_bound = self._bounds(batch_size=self.perm_batch_size, batch=batch)
            # load and transform batch of SNPs
            print("\rCalculate perm test statistics for SNPs %d to %d" % (lower_bound, upper_bound), end='')
            if perm_method == 'y':
                US = self._s_matrix(lower_bound=lower_bound, upper_bound=upper_bound, save_meta=False)  # shape: (n,b)
                # transform data
                US = self.transform_input(X=US, U=self.U)
                # get 3D copy of S for permutations
                US = self.get_3d_copy(v=US, batch_size=self.perm)  # shape: (p,n,b)
            else:
                US = self._s_matrix(lower_bound=lower_bound, upper_bound=upper_bound, device=torch.device("cpu"),
                                    save_meta=False)  # shape: (n,b)
                US = self.permute(data=US)  # shape: (p,n,b)
                # transform data
                US = self.transform_input(X=US, U=self.U)  # shape: (p,n,b)
            # calculate residual sum of squares
            RSS = self.get_rss_perm(S=US, ZD=ZD, ZDZ=ZDZ, v_g=v_g)  # shape: (p,b)
            # calculate test statistics
            test_stat.append(self.get_f_score(rss0=torch.t(RSS_0.repeat(RSS.shape[1], 1)),
                                              rss1=RSS).to(torch.device("cpu")))  # shape: (p,b))
            if self.device.type != "cpu":
                with torch.cuda.device(self.device):
                    del RSS
                    del US
                    torch.cuda.empty_cache()
        test_stat = torch.cat(test_stat, dim=1).to(torch.device("cpu"))  # shape: (p,m)
        time_test_stats = time.time()
        print("\nHave perm test statistics. Elapsed time: ", time_test_stats - var_comp_time)
        if adj_p_value:
            # calculate permutation-based p-values
            self.perm_p_val = self.get_perm_p_value(perm_test_stats=test_stat)  # shape: (m)
            print("Have adjusted p-values")
        # calculate Westfall-Young permutation-based threshold
        self.min_p_value = self.get_min_p_value(test_stat=test_stat)  # shape: (p)
        print("Have minimal p-values. Elapsed time: ", time.time() - time_test_stats)

    def estimate_delta(self, gridlen: int = 100, logdelta_min: int = -10, logdelta_max: int = 10,
                       reml: bool = True) -> torch.tensor:
        """
        Estimate ratio of variance components delta of LMM
        Get grid of evenly divided delta values on logarithmic scale and compute neg loglikelihood for each

        :param gridlen: length of grid, default=100
        :param logdelta_min: lower bound for delta (log value), default=-10
        :param logdelta_max: upper bound for delta (log value), default=10
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: optimal delta
        """
        deltas = torch.exp(torch.linspace(start=logdelta_min, end=logdelta_max, steps=gridlen + 1, device=self.device))
        neglogs = self.negloglikelihood(delta=deltas, Uy=self.Uy, UZ=self.UZ, reml=reml)
        neglogs.to(self.device)
        delta_opt = self._minimize(Uy=self.Uy, UZ=self.UZ, deltas=deltas, neglogs=neglogs, gridlen=gridlen, reml=reml)
        return delta_opt

    def _minimize(self, Uy: torch.tensor, UZ: torch.tensor, deltas: torch.tensor, neglogs: torch.tensor,
                  gridlen: int = 100, reml: bool = True) -> torch.tensor:
        """
        minimize negative loglikelihood function with brent search

        :param Uy: transformed phenotype vector U'y
        :param UZ: transformed vector of fixed effects U'Z
        :param deltas: tensor with possible delta values in ascending order
        :param neglogs: tensor with negative loglikelihood value for each delta
        :param gridlen: length of delta grid, default=100
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: optimal delta
        """
        tmp = torch.argmin(neglogs)
        delta_opt = deltas[tmp]
        neglog_opt = neglogs[tmp]
        # use brent search for each triple in grid
        for i in range(gridlen - 1):
            if (neglogs[i + 1] < neglogs[i]) and (neglogs[i + 1] < neglogs[i + 2]):
                delta_tmp, neglog_tmp, niters = brent.brent_search(f=self.negloglikelihood, a=deltas[i],
                                                                   b=deltas[i + 2], x=deltas[i + 1], fx=neglogs[i + 1],
                                                                   Uy=Uy, UZ=UZ, reml=reml)
                if neglog_tmp < neglog_opt:
                    delta_opt = delta_tmp
                    neglog_opt = neglog_tmp
        return delta_opt

    def negloglikelihood(self, delta: torch.tensor, UZ: torch.tensor, Uy: torch.tensor, reml: bool = True) \
            -> torch.tensor:
        """
        compute negative loglikelihood for one delta value or several values in parallel

        :param delta: ratio of variance components
        :param UZ: transformed fixed effects U'Z
        :param Uy: transformed phenotype U'y
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: negative loglikelihood
        """
        if delta.ndim == 0:
            D = self.D + delta
        else:
            D = self._d_delta(delta=delta, batch_size=len(delta))  # shape: (b,1,n)
        ZD = self._zd(UZ=UZ, D=D)
        ZDZ = self._zdz(UZ=UZ, ZD=ZD)
        beta = self._beta(ZDZ=ZDZ, ZDy=torch.matmul(ZD, Uy))
        sigma = self._sigma(D=D, Uy=Uy, UZ=UZ, beta=beta, reml=reml)
        if D.ndim == 1:
            logdetD = torch.sum(torch.log(D))
        else:
            logdetD = torch.sum(torch.squeeze(torch.log(D)), 1)
        if not reml:
            return (self.dataset.n_samples*torch.log(2*torch.pi*sigma) + logdetD + self.dataset.n_samples) / 2
        else:
            if UZ.ndim == 2:
                logdetZ = torch.logdet(torch.matmul(torch.t(UZ), UZ))
            elif UZ.ndim == 3:
                logdetZ = torch.logdet(torch.matmul(torch.transpose(UZ, dim0=1, dim1=2), UZ))
            else:
                logdetZ = torch.logdet(torch.matmul(torch.transpose(UZ, dim0=2, dim1=3), UZ))
            logdetZDZ = torch.logdet(ZDZ)
            return (self.freedom_deg*torch.log(2*torch.pi*sigma) + logdetD + self.freedom_deg - logdetZ + logdetZDZ) / 2

    def compute_var_components(self, D: torch.tensor, UZ: torch.tensor, ZD: torch.tensor, ZDZ: torch.tensor,
                               reml: bool = True) -> tuple:
        """
        Compute variance components v_g^2 and v_e^2 with Var(y) = v_g^2K + v_e^2I

        :param D: vector with eigenvalues of K
        :param UZ: transformed fixed effects U'Z
        :param ZD: precomputed matrix product of (U'Z)'D^-1
        :param ZDZ: precomputed matrix product of (U'Z)'D^-1(U'Z)
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: v_g^2 and v_e^2
        """
        beta = self._beta(ZDZ=ZDZ, ZDy=torch.matmul(ZD, self.Uy))
        v_g = self._sigma(D=D, Uy=self.Uy, UZ=UZ, beta=beta, reml=reml)
        v_e = self.delta * v_g
        return v_g, v_e

    def get_rss_h0(self, sigma_opt: bool = True, reml: bool = True) -> torch.tensor:
        """
        Compute residual sum of squares of H0 (marker has no effect on phenotype),
        i.e. for fixed effects Z, covariance matrix V and phenotype y compute:
            b = (Z'V^{-1}Z)^{-1}Z'V^{-1}y
            rss = (y-Zb)'V^{-1}(y-Zb)
        note that for optimal sigma_g rss=n-c (REML) or rss=n (ML)

        :param sigma_opt: if True return degrees of freedom, default is True
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: residual sum of squares
        """
        if sigma_opt:
            if reml:
                return torch.tensor(self.dataset.n_samples - self.dataset.fixed.shape[1], device=self.device)
            else:
                return torch.tensor(self.dataset.n_samples, device=self.device)
        else:
            raise NotImplementedError

    def get_rss_and_se(self, D: torch.tensor, S: torch.tensor, ZD: torch.tensor, ZDZ: torch.tensor) -> tuple:
        """
        Compute residual sum of squares of alternative hypothesis (marker has effect on phenotype),
        i.e. for a 3D tensor with batches of fixed effects X and 3D tensor with copies of phenotype y:
            beta = (X'D^{-1}X)^{-1}X'D^{-1}y
            rss = (y-Xbeta)'D^{-1}(y-Xbeta)
        Use block-wise computation for beta, i.e., for computation of beta use the fact that X=[Z,s] for fixed
        effects Z and SNP s.

        :param D: vector with eigenvalues of K + ratio of variance components delta; shape: (n)
        :param S: matrix containing several markers in batches; shape: (n,b)
        :param ZD: precomputed matrix product of (U'Z)'D^-1; shape: (c,n)
        :param ZDZ: precomputed matrix product of (U'Z)'D^-1(U'Z); shape: (c,c)

        :return: residual sum of squares, standard error and effect size in batches
        """
        batch_size = S.shape[1]
        # get (X'D^{-1}X)^{-1}
        SD, XDX = self._xdx(D=D, S=S, ZD=ZD, ZDZ=ZDZ)
        XDX = torch.linalg.pinv(XDX, hermitian=True)
        # compute Z'Dy
        ZDy = self.get_3d_copy(v=torch.matmul(ZD, self.Uy), batch_size=batch_size)  # shape: (b,c,1)
        # compute X'Dy
        SD = torch.matmul(SD, self.Uy).reshape(batch_size, 1, 1)  # shape: (b,1,1)
        # put together 3D tensor
        SD = torch.cat((ZDy, SD), dim=1)  # shape: (b,c+1,1)
        # compute beta
        beta = torch.matmul(XDX, SD)  # shape: (b,c+1,1)
        # compute rss
        S = self._x_batch(X=S, fixed=self.UZ)  # shape (b,n,c+1)
        S = torch.matmul(S, beta)  # shape (b,n,1)
        S = self.get_3d_copy(v=self.Uy, batch_size=batch_size) - S  # shape (b,n,1)
        resD = torch.div(S, torch.unsqueeze(D, 1))
        S = torch.squeeze(torch.matmul(torch.transpose(resD, dim0=1, dim1=2), S)) / self.v_g
        # get standard error
        diag = torch.diagonal(XDX, dim1=1, dim2=2)[:, -1]
        se = torch.sqrt(self.v_g * diag)
        return S, se, torch.squeeze(beta[:, -1])

    def get_f_score(self, rss0: torch.tensor, rss1: torch.tensor) -> torch.tensor:
        """
        Compute tensor of test statistics

        :param rss0: residual sum of squares of H0: marker has no effect on phenotype
        :param rss1: residual sum of squares of H1: marker has effect on phenotype

        :return: F1 score
        """
        return self.freedom_deg * (rss0 - rss1) / rss1

    def get_p_value(self, f_score: float) -> float:
        """
        Compute p-value using survival function of f distribution

        :param f_score: F1 score

        :return: p-value
        """
        return stats.f.sf(f_score, 1, self.freedom_deg)

    # functions for permutations
    def estimate_delta_perm(self, gridlen: int = 100, logdelta_min: int = -10, logdelta_max: int = 10,
                            reml: bool = True) -> torch.tensor:
        """
        Estimate ratio of variance components delta of LMM for permutations
        Get grid of evenly divided delta values on logarithmic scale and compute neg loglikelihood for each

        :param gridlen: length of grid, default=100
        :param logdelta_min: lower bound for delta (log value), default=-10
        :param logdelta_max: upper bound for delta (log value), default=10
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: tensor with optimal delta for each permutation
        """
        deltas = torch.exp(torch.linspace(start=logdelta_min, end=logdelta_max, steps=gridlen + 1, device=self.device))
        if self.UZ.ndim == 2:
            # for perm method y: same U'Z for each permutation
            neglogs = self.negloglikelihood(delta=deltas, Uy=self.get_4d_copy(v=self.Uy, batch_size=len(deltas)),
                                            UZ=self.UZ, reml=reml)
        else:
            # for perm method x: have different U'Z for each permutation
            neglogs = self.negloglikelihood(delta=deltas, Uy=self.get_4d_copy(v=self.Uy, batch_size=len(deltas)),
                                            UZ=self.get_4d_copy(v=self.UZ, batch_size=len(deltas)), reml=reml)
        neglogs.to(self.device)
        delta_opt = []
        if self.UZ.ndim == 2:
            # for perm method y: same U'Z for each permutation
            for i in range(self.perm):
                delta_opt.append(self._minimize(Uy=self.Uy[i, :, 0], UZ=self.UZ, deltas=deltas, neglogs=neglogs[i, :],
                                                gridlen=100, reml=True))
        else:
            # for perm method x: have different U'Z for each permutation
            for i in range(self.perm):
                delta_opt.append(self._minimize(Uy=self.Uy[i, :, 0], UZ=self.UZ[i, :, :], deltas=deltas,
                                                neglogs=neglogs[i, :], gridlen=100, reml=True))
        return torch.tensor(delta_opt, device=self.device)

    def get_rss_perm(self, S: torch.tensor, ZD: torch.tensor, ZDZ: torch.tensor, v_g: torch.tensor) -> torch.tensor:
        """
        Compute residual sum of squares of alternative hypothesis (marker has effect on phenotype) with permutations,
        i.e. for a 4D tensor with copies of batches of fixed effects Z and markers S and 4D tensor with copies of
        permutations of phenotype y:
            b = (X'X)^{-1}X'y
            rss = (y-Xb)'(y-Xb)
        Use block-wise computation for beta, i.e., for computation of beta use the fact that X=[Z,s] for fixed
        effects Z and SNP s.

        :param S: matrix containing batch of markers, shape: (p,n,b)
        :param ZD: 3D tensor containing matrix product (U'Z)'D^{-1} for each permutation, shape: (p,c,n)
        :param ZDZ: 3D tensor containing matrix product (U'Z)'D^{-1}(U'Z) for each permutation, shape: (p,c,c)
        :param v_g: tensor containing genetic variance component for each permutation, shape: (p)

        :return: residual sum of squares in batches
        """
        batch_size = S.shape[2]
        y_batch = self.get_4d_copy(v=self.Uy, batch_size=batch_size)  # shape: (p,b,n,1)
        beta = self._beta_perm(S=S, ZD=ZD, ZDZ=ZDZ, y_batch=y_batch, batch_size=batch_size)  # shape: (p,b,c+1,1)
        # compute residuals
        S = self._x_batch(X=S, fixed=self.UZ)  # shape: (p,b,n,c+1)
        S = y_batch - torch.matmul(S, beta)  # shape: (p,b,n,1)
        # compute residual sum of squares
        rss = torch.div(torch.transpose(S, dim0=2, dim1=3), self.get_4d_copy(v=self.D, batch_size=batch_size))
        rss = torch.squeeze(torch.matmul(rss, S))  # shape: (p,b)
        return torch.t(torch.div(torch.t(rss), torch.unsqueeze(v_g, dim=0)))

    def get_perm_p_value(self, perm_test_stats: torch.tensor) -> torch.tensor:
        """
        Compute permutation-based p-values via
        p = R/(qm) with R being the number of permuted test statistics bigger than the observed test statistic

        :param perm_test_stats: matrix containing test-statistics for all permutations and SNPs, dim (p,m)

        :return: adjusted p-values
        """
        sorted_test_stats, ind = torch.sort(perm_test_stats.flatten())
        n = sorted_test_stats.shape[0]
        test_stats_ind = torch.searchsorted(sorted_test_stats.contiguous(), self.test_stat.contiguous(), right=True)
        adj_p_value = ((n - test_stats_ind) / n).type(torch.float64)
        return torch.where(adj_p_value == 0., 1 / n, adj_p_value)

    def get_min_p_value(self, test_stat: torch.tensor) -> torch.tensor:
        """
        Compute minimal p-values for each permutation:
        First search the maximal test statistic for each permutation, since the survival function is decreasing, this
        gives the minimal p-value

        :param test_stat: matrix containing test-statistics for all permutations and SNPs, dim (p,m)

        :return: vector containing the minimal p-value for each permutation
        """
        max_test_stats, _ = torch.max(test_stat, dim=1)
        min_p_val = []
        for test in max_test_stats:
            min_p_val.append(self.get_p_value(f_score=test))
        return torch.tensor(min_p_val)

    # functions to compute intermediate results
    @staticmethod
    def _zd(UZ: torch.tensor, D: torch.tensor) -> torch.tensor:
        """
        Compute (U'Z)'D^{-1} for fixed effects Z of shape (n,c) or (p,n,c)

        :param UZ: transformed fixed effects U'Z
        :param D: vector with eigenvalues of K + ratio of variance components delta

        :return: Z'D^{-1}
        """
        if UZ.ndim == 2:
            return torch.div(torch.t(UZ), D)
        elif UZ.ndim == 3:
            return torch.div(torch.transpose(UZ, dim0=1, dim1=2), D)
        elif UZ.ndim == 4:
            return torch.div(torch.transpose(UZ, dim0=2, dim1=3), D)

    @staticmethod
    def _zdz(UZ: torch.tensor, ZD: torch.tensor) -> torch.tensor:
        """
        Compute (U'Z)'D^{-1}(U'Z) for fixed effects Z of shape (c,c) or (p,c,c)

        :param UZ: transformed fixed effects U'Z
        :param ZD: precomputed (U'Z)'D^{-1}

        :return: (U'Z)'D^{-1}(U'Z)
        """
        return torch.matmul(ZD, UZ)

    @staticmethod
    def _beta(ZDZ: torch.tensor, ZDy: torch.tensor) -> torch.tensor:
        """
        compute effect size beta = ((U'Z)'D^-1(U'Z))^-1(U'Z)'D^-1(U'y)

        :param ZDZ: precomputed matrix product of (U'Z)'D^-1(U'Z)
        :param ZDy: precomputed matrix product of (U'Z)'D^-1(U'y)

        :return: beta
        """
        return torch.linalg.solve(ZDZ, ZDy)

    def _sigma(self, D: torch.tensor, Uy: torch.tensor, UZ: torch.tensor, beta: torch.tensor, reml: bool = True) \
            -> torch.tensor:
        """
        compute variance component v_g^2 = ((U'y)-(U'Z)beta)'D^-1((U'y)-(U'Z)beta)/(n-c)

        :param D: vector with eigenvalues of K + ratio of variance components delta
        :param Uy: transformed phenotype U'y, shape (n)
        :param UZ: transformed fixed effects U'Z, shape (n,c)
        :param beta: effect size, shape (c)
        :param reml: if True use REML estimate, if False use ML, default=True

        :return: v_g^2
        """
        if D.ndim == 3:
            if Uy.ndim == 1:
                Uy = self.get_3d_copy(v=Uy, batch_size=D.shape[0])
            if beta.ndim == 2:
                beta = torch.unsqueeze(beta, 2)
        res = Uy - torch.matmul(UZ, beta)
        res = torch.multiply(res, res)
        if D.ndim == 1:
            res = torch.sum(torch.div(res, D))
        elif res.ndim == 3:
            res = torch.div(torch.transpose(res, dim0=1, dim1=2), D)
            res = torch.sum(torch.squeeze(res), 1)
        elif res.ndim == 4:
            res = torch.div(torch.transpose(res, dim0=2, dim1=3), D)
            res = torch.sum(torch.squeeze(res), 2)
        if not reml:
            return res / self.dataset.n_samples
        else:
            return res / self.freedom_deg

    def _xdx(self, D: torch.tensor, S: torch.tensor, ZD: torch.tensor, ZDZ: torch.tensor) -> tuple:
        """
        Compute (X'D^{-1}X)^{-1} for X=([Z,s_i],...,[Z,s_{i+b-1}]) of shape (b,n,c+1) for fixed effects Z of shape (n,c)
        and SNPs s_j
        For permutations compute 4D version

        :param D: vector with eigenvalues of K + ratio of variance components delta; shape: (n) or (p,1,n)
        :param S: matrix with batch of b SNPs (n,b) or (p,n,b)
        :param ZD: Z'D^{-1} for fixed effects Z and matrix of eigenvalues+delta D (c,n) or (p,c,n) for perm
        :param ZDZ: Z'D^{-1}Z for fixed effects Z and matrix of eigenvalues+delta D (c,c) or (p,c,c) for perm

        :return: S'D^{-1} and (X'D^{-1}X)^{-1}
        """
        if ZD.ndim == 2:
            batch_size = S.shape[1]
            # compute Z'Ds_i for each SNP s_i in batches
            ZDS = torch.unsqueeze(torch.t(torch.matmul(ZD, S)), dim=2)  # shape: (b,c,1)
            # compute s_iDs_i for all SNPs in batch
            SD = torch.unsqueeze(torch.div(torch.t(S), D), dim=1)  # shape: (b,1,n)
            XDX = torch.bmm(SD, torch.unsqueeze(torch.t(S), dim=2))  # shape: (b,1,1)
            # put together 3D tensor for XDX
            XDX = torch.cat((torch.cat((self.get_3d_copy(v=ZDZ, batch_size=batch_size), ZDS), dim=2),
                             torch.cat((torch.transpose(ZDS, dim0=1, dim1=2), XDX), dim=2)), dim=1)  # shape: (b,c+1,c+1)
        elif ZD.ndim == 3:
            batch_size = S.shape[2]
            # get 4D copy of ZDZ for batch
            ZDZ_4d = self.get_4d_copy(v=ZDZ, batch_size=batch_size)  # shape: (p,b,c,c)
            # compute Z'D^{-1}S
            ZDS = torch.unsqueeze(torch.transpose(torch.matmul(ZD, S), dim0=1, dim1=2), 3)  # shape: (p,b,c,1)
            # compute S'D^{-1}S
            St = torch.transpose(S, dim0=1, dim1=2)  # shape: (p,b,n)
            SD = torch.unsqueeze(torch.divide(St, self.D), dim=2)  # shape: (p,b,1,n)
            # compute S'D^{-1}S
            XDX = torch.matmul(SD, torch.unsqueeze(St, dim=3))
            # put together X'D^{-1}X
            XDX = torch.concat((torch.transpose(ZDS, dim0=2, dim1=3), XDX), dim=3)
            XDX = torch.concat((torch.concat((ZDZ_4d, ZDS), dim=3), XDX), dim=2)  # shape: (p,b,c+1,c+1)
        else:
            raise Exception('Can only compute XDX for 2D or 3D version of ZD.')
        return SD, XDX

    def _beta_perm(self, S: torch.tensor, ZD: torch.tensor, ZDZ: torch.tensor, y_batch: torch.tensor, batch_size: int) \
            -> torch.tensor:
        """
        Compute betas for permutations in 4D tensor using block-wise computations

        :param S: matrix containing batch of markers, shape: (p,n,b)
        :param ZD: 3D tensor containing matrix product (U'Z)'D^{-1} for each permutation, shape: (p,c,n)
        :param ZDZ: 3D tensor containing matrix product (U'Z)'D^{-1}(U'Z) for each permutation, shape: (p,c,c)
        :param y_batch: 4D copy of permutations of phenotype vector, shape: (p,b,n,1)
        :param batch_size: number of markers

        :return: 4D tensor with beta values for all markers nad permutations, shape: (p,b,c+1,1)
        """
        # get S'D^{-1}S and X'D^{-1}X
        SD, XDX = self._xdx(D=self.D, S=S, ZD=ZD, ZDZ=ZDZ)  # shape: (p,b,1,n), (p,b,c+1,c+1)
        # get X'D^{-1}y
        XDy = self.get_4d_copy(v=torch.matmul(ZD, self.Uy), batch_size=batch_size)  # shape: (p,b,c,1)
        SD = torch.matmul(SD, y_batch)  # shape: (p,b,1,1)
        XDy = torch.concat((XDy, SD), dim=2)  # shape: (p,b,c+1,1)
        # get beta of shape: (p,b,c+1,1)
        return self._beta(ZDZ=XDX, ZDy=XDy)

    # functions for data transformation
    @staticmethod
    def transform_input(X: torch.tensor, U: torch.tensor) -> torch.tensor:
        """
        compute U'X

        :param X: input vector/matrix
        :param U: input matrix

        :return: product with transpose
        """
        return torch.matmul(torch.t(U), X)

    def _d_delta(self, delta: torch.tensor, batch_size: int):
        """
        get 3D tensor with D + delta*I as batches for diagonal matrix with eigenvalues D and different variance
        component ratios delta. If delta is one value, return tensor with b copies of D+delta.

        :param delta: variance component ratio shape: (b) or (1)
        :param batch_size: number of needed copies of D

        :return: D + delta of shape (b,1,n)
        """
        if delta.ndim == 1:
            return torch.unsqueeze(self.D.repeat(batch_size, 1) + torch.unsqueeze(delta, 1), 1)
        else:
            return torch.unsqueeze((self.D + delta).repeat(batch_size, 1), 1)

    def _s_matrix(self, lower_bound: int, upper_bound: int, device=None, save_meta: bool = True) -> torch.tensor:
        """
        load batch of markers to specified device

        :param lower_bound: lower bound of marker batch
        :param upper_bound: upper bound of marker batch
        :param device: either cpu or cuda device
        :param save_meta: if genotype is loaded batch-wise, set to False for permutations to prevent saving of meta info

        :return: matrix with markers of shape (n,upper_bound-lower_bound)
        """
        if device is None:
            device = self.device
        if self.dataset.X is None:
            # load X batch-wise
            self.dataset.load_genotype_batch_wise(device=device, save_meta=save_meta, snp_lower_index=lower_bound,
                                                  snp_upper_index=upper_bound)  # shape: (n,b)
            S = self.dataset.X  # shape: (n,b)
            self.dataset.reset_genotype()
        else:
            # get X_batch if X was completely loaded before
            S = self.dataset.X[:, lower_bound:upper_bound].to(device)  # shape: (n,b)
        return S

    def _x_batch(self, X: torch.tensor, fixed: torch.tensor) -> torch.tensor:
        """
        Create 3D or 4D tensor where each matrix in the 3D tensor contains the same fixed effects and a different SNP,
        and the 4D tensor contains copies of the 3D tensors

        :param X: genotype matrix/tensor of shape (n,b) or (p,n,b)
        :param fixed: matrix/tensor of fixed effects of shape (n,c) or (p,n,c)

        :return: tensor of shape (b,n,c+1) or (p,b,n,c+1)
        """
        if X.ndim == 2:
            b = self.get_3d_copy(v=fixed, batch_size=X.shape[1])
            return torch.cat((b, torch.transpose(torch.unsqueeze(X, 0), 0, 2)), dim=2)
        elif X.ndim == 3:
            b = self.get_4d_copy(v=fixed, batch_size=X.shape[2])
            return torch.cat((b, torch.unsqueeze(torch.transpose(X, dim0=1, dim1=2), 3)), dim=3)

    @staticmethod
    def get_3d_copy(v: torch.tensor, batch_size: int) -> torch.tensor:
        """
        Create 3D tensor with copies of input tensor

        :param v: vector/matrix of shape (n) or (n,c)
        :param batch_size: batch size of new 3D tensor

        :return: tensor of copies of v with shape (batch_size,n,1) or (batch_size,n,c)
        """
        if v.ndim == 1:
            return torch.unsqueeze(v.expand(batch_size, v.shape[0]), 2)
        if v.ndim == 2:
            return v.expand(batch_size, v.shape[0], v.shape[1])

    @staticmethod
    def get_4d_copy(v: torch.tensor, batch_size: int) -> torch.tensor:
        """
        Create 4D tensor with copies of input tensor

        :param v: tensor of shape (p,n,c)
        :param batch_size: batch size of new 4D tensor

        :return: tensor of copies of v with shape (p,b,n,c)
        """
        return torch.transpose(v.expand(batch_size, v.shape[0], v.shape[1], v.shape[2]), dim0=0, dim1=1)

    # helper functions
    def _bounds(self, batch_size: int, batch: int) -> tuple:
        """
        compute upper and lower bound for natch-wise computations

        :param batch_size: number of markers within batch
        :param batch: number of batch

        :return: lower and upper bound
        """
        lower_bound = batch * batch_size
        upper_bound = (batch + 1) * batch_size
        if upper_bound > self.dataset.n_snps:
            upper_bound = self.dataset.n_snps
        return lower_bound, upper_bound
