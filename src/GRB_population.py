import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import yaml
from pathlib import Path
from io_grb_pop import read_column, root_dir
from cosmology import init_cosmology, Lum_dist
from observational_constraints import load_observational_constraints
import observational_constraints as obs
import stats as st
import physics as ph

log = logging.getLogger(__name__)

try:
    import f90_functions as f90f
    f90 = True
except ImportError:
    f90 = False
    log.error("Could not import f90_functions, f90 set to False")


class GRBPopulation:
    def __init__(self, Nb_GRBs, output_dir=None):
        self.Nb_GRBs = Nb_GRBs
        self.properties = pd.DataFrame({})
        self.mock_constraints = {}
        self.likelihood_params = {}
        if output_dir is None:
            self.output_dir = Path().cwd()
        else:
            self.output_dir = output_dir

    def draw_L(self, Nb_GRBs=None, model='EPL', z=None, run_mode=None, savefig=False, **params):
        """ Draw L from a power law """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs
        if model == 'EPL':
            logLmin = params['logLmin']
            logLmax = params['logLmax']
            slope = params['slope']
            k_evol = params['k_evol']
            t = np.random.rand(Nb_GRBs)
            L = 10.**logLmin*(1. - (1. - (10.**(logLmax-logLmin))**(1.-slope))*t)**(1. / (1.-slope))
            if k_evol != 0.0:
                if z is None:
                    try:
                        z = self.properties['z']
                        L *= (1 + z)**k_evol
                    except KeyError:
                        raise ValueError('Could not find z to use in L_draw')
                else:
                    L *= (1 + z)**k_evol
        elif model == 'ES':
            logLmin = params['logLmin']
            logLbreak = params['logLbreak']
            logLmax = logLbreak + 1.5  # to avoid extremely low probabilities beyond logLbreak + 1.5
            slope = params['slope']
            k_evol = params['k_evol']
            logL_range = np.linspace(logLmin, logLmax, 1000)
            logL_pdf_unnormed = self.Schechter_log(logL_range, logLbreak=logLbreak, slope=slope)
            logL_pdf = logL_pdf_unnormed / logL_pdf_unnormed.sum()

            try:
                L = 10.**(np.random.choice(logL_range, size=Nb_GRBs, p=logL_pdf))
            except Exception as e:
                L = np.zeros(Nb_GRBs)
                log.warning(f'Raised exception {e}')

            if k_evol != 0.0:
                if z is None:
                    if 'z' not in self.properties.keys():
                        raise KeyError('Could not find z to use in draw_L')
                    z = self.properties['z']
                    L *= (1 + z)**k_evol
                else:
                    L *= (1 + z)**k_evol
        elif model == 'EBPL':
            raise NotImplementedError
        else:
            raise ValueError("Invalid model for L drawing")

        self.properties['L'] = L

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting L pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            ax.hist(np.log10(L),
                    bins=np.arange(logLmin, logLmax, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(L [erg/s])')
            ax.set_ylabel('L PDF')
            ax.set_yscale('log')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'L_draw.pdf')

        return L.copy()

    def draw_z(self, cosmo=None, Nb_GRBs=None, zmax=20, model='SH', run_mode=None, savefig=False, **params):
        """ Draw z from a redshift distribution """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if cosmo is None:
            try:
                cosmo = init_cosmology(Path().cwd().parent/'data'/'cosmology')
            except FileNotFoundError:
                raise ValueError("Could not draw z because no cosmology provided")

        redshift = cosmo['redshift']
        dVdz = cosmo['dVdz']

        if zmax > redshift[-1]:
            log.warning("In draw_z, zmax is greater than largest value in redshift table")
        z_range = redshift[redshift <= zmax]
        comoving_volume = dVdz[redshift <= zmax] / (1. + z_range)

        if model == 'SH':
            zm = params['zm']
            z_a = params['a']
            z_b = params['b']
            z_pdf_comov = self.SH(z_range, a=z_a, b=z_b, zm=zm)
            z_pdf_unnormed = z_pdf_comov * comoving_volume
        elif model == 'BExp':
            zm = params['zm']
            z_a = params['a']
            z_b = params['b']
            z_pdf_comov = self.GRBrate_exp(z_range, a=z_a, b=z_b, zm=zm)
            z_pdf_unnormed = z_pdf_comov * comoving_volume
        else:
            raise NotImplementedError

        z_pdf = z_pdf_unnormed / z_pdf_unnormed.sum()

        try:
            z = np.random.choice(z_range, size=Nb_GRBs, p=z_pdf)
        except Exception as e:
            z = np.zeros(Nb_GRBs)
            log.warning(f'Raised exception {e}')

        self.properties['z'] = z

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting z pdf")
            fig, axes = plt.subplots(2, figsize=(6, 7), tight_layout=True)
            delta_z = z_range[1]-z_range[0]
            axes[0].hist(z, bins=50, density=True, color='lightgray', label='Drawings',
                         edgecolor='k', linewidth=0.5)
            axes[0].plot(z_range, z_pdf/delta_z, label='PDF')
            axes[0].set_ylabel('z PDF')
            axes[0].set_yscale('log')
            axes[0].legend()
            axes[1].plot(z_range, z_pdf_comov, label='comoving rate')
            axes[1].set_ylabel(r'$\rm[yr{-1}\,Mpc^{-3}]$')
            axes[1].legend()
            axes[1].set_xlabel('Redshift')
            axes[1].set_yscale('log')
            if savefig:
                self.save_fig(fig, 'z_draw.pdf')
        return z.copy()

    def draw_Ep(self, Nb_GRBs=None, model='LN', L=None, run_mode=None, savefig=False, **params):
        """ Draw Ep from LogNormal distribution """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs
        Ep0 = params["Ep0"]
        sigmaEp = params["sigmaEp"]
        alphaA = params["alpha_amati"]
        L0 = 1.6e52  # erg/s

        if model == 'LN':
            t = np.random.normal(np.log10(Ep0), sigmaEp, Nb_GRBs)
            Ep = 10.**t
        elif model == 'A':
            if L is None:
                try:
                    L = self.properties['L']
                except KeyError:
                    raise ValueError("Could not find L for Amati Ep drawing")
            t = np.random.normal(0, sigmaEp, Nb_GRBs)
            Ep = Ep0 * (L/L0)**alphaA * 10.**(t * np.sqrt(1. + alphaA**2))
        else:
            raise ValueError("Invalid model for Ep drawing")

        self.properties['Ep'] = Ep

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting Ep pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            logEpmin = np.log10(Ep0) - 5*sigmaEp
            logEpmax = np.log10(Ep0) + 5*sigmaEp
            ax.hist(np.log10(Ep),
                    bins=np.arange(logEpmin, logEpmax, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(Ep [keV])')
            ax.set_ylabel('Ep PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'Ep_draw.pdf')
        return Ep.copy()

    def draw_spec(self, Nb_GRBs=None, model='FBand', run_mode=None, savefig=False, data_dir=root_dir/'data', **params):
        """
            Draw the spectral parameters alpha, beta and ktild
        """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if model == 'FBand':
            _ktild = f90f.f90f.calc_ktild(alpha=params['alpha'], beta=params['beta'])
            ktild = _ktild * np.ones(Nb_GRBs)
            alpha = params['alpha'] * np.ones(Nb_GRBs)
            beta = params['beta'] * np.ones(Nb_GRBs)
        elif (model == 'GBM_Band_old') or (model == 'GBM_Band'):
            if (model == 'GBM_Band_old'):
                alpha_file = 'alpha_GBM.txt'
                beta_file = 'beta_GBM.txt'
            elif (model == 'GBM_Band'):
                alpha_file = 'good_alpha_GBM.txt'
                beta_file = 'good_beta_GBM.txt'
            alpha = self.draw_from_cdf_file(data_dir/alpha_file, N_draws=Nb_GRBs)
            beta = self.draw_from_cdf_file(data_dir/beta_file, N_draws=Nb_GRBs)
            # alpha_range, alpha_pdf_unnormed = self.create_pdf_from_cdf(data_dir/alpha_file)
            # beta_range, beta_pdf_unnormed = self.create_pdf_from_cdf(data_dir/beta_file)
            # alpha_pdf = alpha_pdf_unnormed/alpha_pdf_unnormed.sum()
            # beta_pdf = beta_pdf_unnormed/beta_pdf_unnormed.sum()
            # alpha = np.random.choice(alpha_range, size=Nb_GRBs, p=alpha_pdf)
            # beta = np.random.choice(beta_range, size=Nb_GRBs, p=beta_pdf)

            if (model == 'GBM_Band_old'):  # because of wrong convention
                alpha = -alpha
                beta = -beta
                beta = np.where((beta == 2), 2.01, beta)
            ktild = f90f.f90f.calc_ktild(alpha=alpha, beta=beta)
        else:
            raise NotImplementedError

        self.properties['alpha'] = alpha
        self.properties['beta'] = beta
        self.properties['ktild'] = ktild

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting spectral params pdf")
            fig, axar = plt.subplots(2, figsize=(6, 7), tight_layout=True)
            for i, item, name in zip([0, 1], [alpha, beta], ['alpha', 'beta']):
                axar[i].hist(item,
                             bins=30,
                             density=True, color='lightgray',
                             label='Drawings', edgecolor='k', linewidth=0.5)
                axar[i].set_xlabel(name)
                axar[i].set_ylabel(name+' PDF')
                axar[i].legend()
            if savefig:
                self.save_fig(fig, 'alpha_beta_draw.pdf')
        return alpha.copy(), beta.copy(), ktild.copy()

    def draw_t90obs(self, Nb_GRBs=None, run_mode=None, savefig=False, **params):
        """ Draw t90obs from LogNormal distribution """
        log.warning("Jesse, you are drawing from T90obs instead of T90intr."
                    " Make sure you know what you're doing.")
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs
        mu_t90obs = params['mu']
        sigma_t90obs = params['sigma']
        t90obs = 10.**(np.random.normal(mu_t90obs, sigma_t90obs, Nb_GRBs))

        self.properties['t90obs'] = t90obs

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting t90obs pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            logt90obsmin = mu_t90obs - 5*sigma_t90obs
            logt90obsmax = mu_t90obs + 5*sigma_t90obs
            ax.hist(np.log10(t90obs),
                    bins=np.arange(logt90obsmin, logt90obsmax, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(t90obs [s])')
            ax.set_ylabel('t90obs PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 't90obs_draw.pdf')
        return t90obs.copy()

    def draw_t90(self, z_med=None, Nb_GRBs=None, run_mode=None, savefig=False, **params):
        """
            Draw t90 from LogNormal distribution corrected by a median
            redshift
        """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs
        if z_med is None:
            try:
                z_med = params['z_med']
            except KeyError:
                z_med = 1.7
                log.warning("In draw_t90, no z_med was found to correct the t90 distribution."
                            " A default value of z_med = 1.7 was used but you should check this.")
        # Correct the observed mu with the median redshift of the
        # sample on which it was derived
        mu_t90 = params['mu']-np.log10(1.+z_med)
        sigma_t90obs = params['sigma']
        t90 = 10.**(np.random.normal(mu_t90, sigma_t90obs, Nb_GRBs))

        self.properties['t90'] = t90

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting t90 pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            logt90min = mu_t90 - 5*sigma_t90obs
            logt90max = mu_t90 + 5*sigma_t90obs
            ax.hist(np.log10(t90),
                    bins=np.arange(logt90min, logt90max, 0.1),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(t90 [s])')
            ax.set_ylabel('t90 PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 't90_draw.pdf')
        return t90.copy()

    def draw_Cvar(self, Nb_GRBs=None, t90obs=None, run_mode=None, savefig=False, **params):
        """ Draw Cvar from t90obs correlated distribution """
        if Nb_GRBs is None:
            Nb_GRBs = self.Nb_GRBs

        if t90obs is None:
            self._check_properties(necessary_prop=['t90obs'],
                                   func_name='Cvar')
            t90obs = self.properties['t90obs']

        mu_Cvar = params['mu']
        sigma_Cvar = params['sigma']
        correl_slope = params['correl_slope']
        t = np.random.normal(mu_Cvar, sigma_Cvar, Nb_GRBs)
        Cvar = 10.**(t - correl_slope * np.log10(t90obs))
        Cvar = np.where((Cvar >= 1), np.ones(Cvar.shape), Cvar)

        self.properties['Cvar'] = Cvar

        if run_mode == 'debug':
            log.info("Debug mode activated; plotting Cvar pdf")
            fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
            ax.hist(np.log10(Cvar),
                    bins=np.arange(-2, 0.1, 0.05),
                    density=True, color='lightgray',
                    label='Drawings', edgecolor='k', linewidth=0.5)
            ax.set_xlabel('log(Cvar)')
            ax.set_ylabel('Cvar PDF')
            ax.legend()
            if savefig:
                self.save_fig(fig, 'Cvar_draw.pdf')
        return Cvar.copy()

    def draw_GRB_properties(self, cosmo, params, run_mode=None, savefig=False):
        """
            Draw the various GRB properties used in the GRB population
            code. Also calculate some others from the drawn properties.
        """

        self.draw_z(**params['redshift_distribution'], cosmo=cosmo, run_mode=run_mode, savefig=savefig)
        self.draw_L(**params['luminosity_function'], run_mode=run_mode, savefig=savefig)
        self.draw_Ep(**params['peak_energy_distribution'], run_mode=run_mode, savefig=savefig)
        self.draw_spec(**params['spectral_shape'], run_mode=run_mode, savefig=savefig)
        D_L = Lum_dist(self.properties['z'], cosmo)
        Epobs = self.properties['Ep']/(1. + self.properties['z'])

        # The median redshift of the population of GBM_bright is needed
        # to calculate t90 (which is defined as pht_pflx in BATSE band
        # above 0.9 ph/s/cm2)
        self.calc_peak_photon_flux(incl_instruments=['BATSE'])
        GBM_bright_selection = self.properties['pht_pflx_BATSE'] >= 0.9
        z_med = np.median(self.properties[GBM_bright_selection]['z'])
        self.draw_t90(**params['t90obs_distribution'], z_med=z_med, run_mode=run_mode, savefig=savefig)
        t90obs = self.properties['t90'] * (1. + self.properties['z'])
        self.draw_Cvar(**params['Cvar_distribution'], run_mode=run_mode, savefig=savefig)

        Eiso = self.properties['L'] * self.properties['Cvar'] * self.properties['t90']

        self.properties['D_L'] = D_L
        self.properties['Epobs'] = Epobs
        self.properties['Eiso'] = Eiso
        self.properties['t90obs'] = t90obs

        if run_mode == 'debug':
            summary = self.summary()
            print(summary)

        return self.properties

    def draw_GRB_properties_for_MCMC(self, cosmo, params, incl_instruments):
        """
            Draw only the core properties necessary to run the MCMC
            exploration
        """

        self.draw_z(**params['redshift_distribution'], cosmo=cosmo)
        self.draw_L(**params['luminosity_function'])
        self.draw_Ep(**params['peak_energy_distribution'])
        self.draw_spec(**params['spectral_shape'])

        D_L = Lum_dist(self.properties['z'], cosmo)
        Epobs = self.properties['Ep']/(1. + self.properties['z'])
        self.properties['D_L'] = D_L
        self.properties['Epobs'] = Epobs
        # Use the function from physics module to avoid unecessary
        # checks present in the class method
        ph.calc_peak_photon_flux(GRB_prop=self.properties,
                                 incl_instruments=incl_instruments)

        return self.properties

    def draw_from_cdf_file(self, filename, N_draws, **args):
        """
            Draw from an ascii file that contains two columns:
            x, CDF(x)
        """
        value_range = read_column(filename, 0)
        cdf = read_column(filename, 1)
        draws = np.random.rand(N_draws)
        values = value_range[cdf.searchsorted(draws)]
        return values

    def calc_peak_photon_flux(self, incl_instruments):
        """
            Calculate the peak photon flux for every GRB in the
            population.
        """
        necessary_prop = ['L', 'z', 'Ep', 'alpha', 'beta', 'D_L']
        self._check_properties(necessary_prop, func_name='peak photon flux')
        ph.calc_peak_photon_flux(GRB_prop=self.properties,
                                 incl_instruments=incl_instruments)
        return

    def calc_peak_energy_flux(self, incl_instruments):
        """
            Calculate the peak energy flux for every GRB in the
            population.
        """
        necessary_prop = ['L', 'z', 'Ep', 'alpha', 'beta', 'D_L']
        self._check_properties(necessary_prop, func_name='peak energy flux')
        ph.calc_peak_energy_flux(GRB_prop=self.properties,
                                 incl_instruments=incl_instruments)
        return

    def calc_photon_fluence(self, incl_instruments):
        """
            Calculate the photon fluence in units of ph/cm2 over the T90
            of the burst.
        """
        necessary_prop = ['t90obs', 'Cvar']
        self._check_properties(necessary_prop, func_name='photon fluence')
        ph.calc_photon_fluence(GRB_prop=self.properties,
                               incl_instruments=incl_instruments)
        return

    def calc_energy_fluence(self, incl_instruments):
        """
            Calculate the photon fluence in units of ph/cm2 over the T90
            of the burst.
        """
        necessary_prop = ['t90obs', 'Cvar']
        self._check_properties(necessary_prop, func_name='energy fluence')
        ph.calc_energy_fluence(GRB_prop=self.properties,
                               incl_instruments=incl_instruments)
        return

    def calc_det_prob(self, incl_samples, **ECLAIRs_args):
        """
            Calculates the detection probability for the included
            samples
        """
        ph.calc_det_prob(GRB_prop=self.properties,
                         incl_samples=incl_samples, **ECLAIRs_args)
        return

    def _check_properties(self, necessary_prop, func_name):
        """
            Convenience, private function to check if the inputs
            required for a given calculation are present in the
            population properties
        """
        for prop in necessary_prop:
            if prop not in self.properties.keys():
                raise KeyError(f"Property {prop} must exist before you attempt to calculate the"
                               f" {func_name}.")

    def create_pdf_from_cdf(self, filename, **args):
        x = read_column(filename, 0, **args)
        cdf = read_column(filename, 1, **args)
        pdf = np.zeros(len(cdf))
        pdf = cdf[1:]-cdf[:-1]
        for i in range(1,len(cdf)):
            pdf[i] = cdf[i]-cdf[i-1]
        return x, pdf

    def Schechter_log(self, logL, logLbreak, slope):
        """ Returns the unnormalized Schechter function
            Expects Lum arguments to be in log scale """
        x = 10.**(logL - logLbreak)
        Sch = x**(1.-slope) * np.exp(-x)
        return Sch

    def SH(self, z, zm=2, a=2.37, b=1.8, nu=0.178, IMF_norm=0.007422):
        """
            Springel-Hernquist+03 functional form for the cosmic SFR.
            Default are given the values of Vangioni+15.
            Returns an event rate in units of yr-1 Mpc-3
            Note : nu is in units of Msun/yr/Mpc3 and IMF_norm in units of Msun-1
        """
        return IMF_norm * nu * a * np.exp(b*(z-zm)) / ((a-b) + b*np.exp(a*(z-zm)))

    def GRBrate_exp(self, z, a=1.1, b=-0.57, zm=1.9, norm=0.00033313):
        """
            GRB rate as parametrized by a broken exponential function.
            Default values are chosen as best fit from SFR of Vangioni+15
            Normalization is done on the same SFR, yielding units of yr-1 Mpc-3
        """
        if isinstance(z, np.ndarray):
            rate = np.where(z <= zm,
                            np.exp(a*z),
                            np.exp(b*z) * np.exp((a-b)*zm))
        else:
            if z <= zm:
                rate = np.exp(a*z)
            else:
                rate = np.exp(b*z) * np.exp((a-b)*zm)
        return norm*rate

    def create_mock_constraint(self, obs_constraints=None):
        """
            Create the mock constraints from the current population.
            obs_constraints is a dictionary with the necessary
            information about each constraint.
        """

        if obs_constraints is None:
            with open(root_dir/'init/obs_constraints.yml', 'r') as f:
                obs_constraints = yaml.safe_load(f)
            load_observational_constraints(obs_constraints)

        for name, constraint in obs_constraints.items():
            bins = constraint['bins']
            val_min = constraint['val_min']
            prop_min = constraint['prop_min']
            quantity = constraint['quantity']
            cond = self.properties[prop_min] >= val_min
            mod, _u = np.histogram(self.properties[cond][quantity], bins=bins)
            self.mock_constraints[name] = {'val_min':val_min,
                                           'prop_min':prop_min,
                                           'quantity':quantity,
                                           'bins':bins,
                                           'hist_unnormed':mod,
                                           'err':np.sqrt(mod)}
        return

    def compare_to_observational_constraints(self, obs_constraints, method='chi2'):
        """
            First normalize the population to the observational
            constraints then calculate the likelihood of the current
            population.
        """
        lnL_tot = 0.0
        chi2_tot = 0.0
        for name, constraint in obs_constraints.items():
            # Normalize to observations
            model = self.mock_constraints[name]['hist_unnormed']
            norm = obs.normalize_to_constraint(mod=model,
                                               obs=constraint['hist'],
                                               err=constraint['err'])
            self.mock_constraints[name]['norm'] = norm
            self.mock_constraints[name]['hist'] = norm * model
            # Calculate likelihood
            if method == 'chi2':
                chi2 = st.chi2(mod=self.mock_constraints[name]['hist'],
                               obs=constraint['hist'],
                               err=constraint['err'])
                self.likelihood_params['_'.join(['chi2',name])] = chi2
                chi2_tot += chi2
                lnL = -0.5 * chi2
            elif method == 'pBIL':
                lnL = st.pBIL(mod=self.mock_constraints[name]['hist'],
                              obs=constraint['hist'],
                              sum_ln_oi_factorial=constraint['sum_ln_oi_factorial'])
                # If using pBIL, add a factor of 10 weight to the eBAT6
                # constraint to make it more impactful
                # if name == 'eBAT6':
                #     lnL *= 10
            self.likelihood_params['_'.join(['lnL',name])] = lnL
            lnL_tot += lnL

        self.likelihood_params['lnL_tot'] = lnL_tot
        if method == 'chi2':
            self.likelihood_params['chi2_tot'] = chi2_tot

        return

    def save_all_GRBs(self, output_dir):
        """
            Save all GRBs to a data frame format openable with pandas.
            This file might be big depending on the number of GRBs.
        """
        df = pd.DataFrame(self.properties)
        outfile = output_dir/"GRB_Properties"
        df.to_msgpack(outfile)

        log.info(f"Saved all GRB properties in {outfile}")

        return

    def save_fig(self, fig, filename):
        try:
            fig.savefig(self.output_dir/filename)
        except FileNotFoundError:
            self.output_dir.mkdir()
            log.info(f"Created '{self.output_dir}' to save {filename} ")
            fig.savefig(self.output_dir/filename)
        return

    def summary(self, full_width=80, cell_width=None, sides='|', thick_line="=", thin_line='-'):
        """
            A short summary of the current LGRB population returned as
            a string.
        """
        if full_width % 2 != 0:
            full_width += 1
            log.warning("In GRBPopulation.summary(): the width you are asking for is odd."
                        " It was increased by 1 for esthetic purposes.")
        max_len = max(map(len, self.properties)) + 1
        if cell_width is None:
            cell_width = np.max([max_len, 14])
        center_width = full_width - 2*len(sides)
        half_width = int(center_width / 2)
        summary = "\n" + full_width * thick_line + "\n"
        summary += sides + "SUMMARY".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        # General configuration
        _Nb_GRBs_string = f"Nb_GRBs =".rjust(half_width)
        _Nb_GRBs = f" {self.Nb_GRBs:.2e}".ljust(half_width+1)
        _output_dir_string = f"Output directory =".rjust(half_width)
        _output_dir = f" {self.output_dir.stem}".ljust(half_width+1)
        summary += sides + _Nb_GRBs_string + _Nb_GRBs + sides + "\n"
        summary += sides + _output_dir_string + _output_dir + sides + "\n"
        summary += full_width * thick_line + "\n"
        # Properties
        summary += sides + "Properties".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        summary += '|'.join([' prop'.ljust(cell_width),
                             'median'.center(cell_width),
                             'stdev'.center(cell_width),
                             'min'.center(cell_width),
                             'max'.center(cell_width)]) + "\n"
        summary += full_width * thin_line + "\n"
        for key, prop in self.properties.items():
            med = np.median(prop)
            std = np.std(prop)
            _min = prop.min()
            _max = prop.max()
            _summary = '|'.join([f" {key}".ljust(cell_width),
                                 f"{med:< .4e}".center(cell_width),
                                 f"{std:< .4e}".center(cell_width),
                                 f"{_min:< .4e}".center(cell_width),
                                 f"{_max:< .4e}".center(cell_width)])
            summary += _summary.ljust(full_width) + "\n"

        summary += full_width * thick_line + "\n"
        summary += sides + "Likelihood".center(center_width) + sides + "\n"
        summary += full_width * thick_line + "\n"
        # Likelihood
        for key, val in sorted(self.likelihood_params.items()):
            _key = f"{key} =".rjust(half_width)
            _val = f" {val:.5f}".ljust(half_width+1)
            _summary = sides + _key + _val + sides + "\n"
            summary += _summary

        summary += full_width * thick_line + "\n"

        return summary
