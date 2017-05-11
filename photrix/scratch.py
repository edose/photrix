import numpy as np
import pandas as pd
import statsmodels.regression.mixed_linear_model as sm  # statsmodels version >= 0.8 !

import matplotlib.pyplot as plt

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def argh4():
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gsp

    x1 = np.linspace(0.0, 5.0)
    x2 = np.linspace(0.0, 2.0)

    y1 = np.cos(2 * np.pi * x1) * np.exp(-1 * x1)
    y2 = np.cos(2 * np.pi * x2)

    fig = plt.figure()  # (width, height) in inches
    fig.set_size_inches([25, 5], forward=True)
    # gs = gsp.GridSpec(nrows=1, ncols=4)

    # plt.subplot(4, 1, 1)
    # plt.plot(x1, y1, 'ko-')
    # plt.title('A tale of 2 subplots')
    # plt.ylabel('Damped oscillation')
    plot1(x1, y1)

    plt.subplot(1, 4, 2)
    plt.plot(x2, y2, 'r.-')
    plt.xlabel('time (s)')
    plt.ylabel('Undamped')

    plt.subplot(1, 4, 3)
    plt.plot(x1, y1, 'ko-')
    plt.title('A tale of 2 subplots')
    plt.ylabel('Damped oscillation')

    plt.subplot(1, 4, 4)
    plt.plot(x2, y2, 'r.-')
    plt.xlabel('time (s)')
    plt.ylabel('Undamped')

    plt.show()

    fig = plt.figure()  # (width, height) in inches
    fig.set_size_inches([25, 5], forward=True)
    # gs = gsp.GridSpec(nrows=1, ncols=4)

    # plt.subplot(4, 1, 1)
    # plt.plot(x1, y1, 'ko-')
    # plt.title('A tale of 2 subplots')
    # plt.ylabel('Damped oscillation')
    plot1(x1, y1)
    plt.show()


def plot1(x1, y1):
    plt.subplot(1, 4, 1)
    plt.plot(x1, y1, 'ko-')
    plt.title('A tale of 2 subplots')
    plt.ylabel('Damped oscillation')


def alt():
    import numpy as np
    import matplotlib.pyplot as plt

    # generate some data
    x = np.arange(0, 10, 0.2)
    y = np.sin(x)

    # plot it
    fig = plt.figure(figsize=(4, 16))  # (width, height)
    f, (a0, a1, a2, a3) = plt.subplots(4, 1)  # will need aspect ratio keyword here
    a0.plot(x, y)
    a1.plot(y, x)
    a2.plot(x, y)
    a3.plot(y, x)

    # f.tight_layout()
    plt.show()


def widget():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)
    t = np.arange(0.0, 1.0, 0.001)
    a0 = 5
    f0 = 3
    s = a0 * np.sin(2 * np.pi * f0 * t)
    l, = plt.plot(t, s, lw=2, color='red')
    plt.axis([0, 1, -10, 10])

    axcolor = 'lightgoldenrodyellow'
    axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    axamp = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

    sfreq = Slider(axfreq, 'Freq', 0.1, 30.0, valinit=f0)
    samp = Slider(axamp, 'Amp', 0.1, 10.0, valinit=a0)

    def update(val):
        amp = samp.val
        freq = sfreq.val
        l.set_ydata(amp * np.sin(2 * np.pi * freq * t))
        fig.canvas.draw_idle()

    sfreq.on_changed(update)
    samp.on_changed(update)

    resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

    def reset(event):
        sfreq.reset()
        samp.reset()

    button.on_clicked(reset)

    rax = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
    radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)

    def colorfunc(label):
        l.set_color(label)
        fig.canvas.draw_idle()

    radio.on_clicked(colorfunc)

    plt.show()


def w2():
    from matplotlib.widgets import Cursor
    import numpy as np
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, facecolor='#FFFFCC')

    x, y = 4 * (np.random.rand(2, 100) - .5)
    ax.plot(x, y, 'o')
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)

    # set useblit = True on gtkagg for enhanced performance
    # cursor = Cursor(ax, useblit=True, color='red', linewidth=2)
    cursor = Cursor(ax, color='red', linewidth=2)  # ???

    plt.show()


def toggle():
    # ---> (EVD) I think we could adapt this to using arrow keys to switch between make_plots.

    """ toggle between two images by pressing "t"

    The basic idea is to load two images (they can be different shapes) and plot
    them to the same axes.  Then, toggle the visible property of
    them using keypress event handling

    If you want two images with different shapes to be plotted with the same
    extent, they must have the same "extent" property

    As usual, we'll define some random images for demo.  Real data is much more
    exciting!

    Note, on the wx backend on some platforms (e.g., linux), you have to
    first click on the figure before the keypress events are activated.
    If you know how to fix this, please email us!

    """

    import matplotlib.pyplot as plt
    import numpy as np

    # two images x1 is initially visible, x2 is not
    x1 = np.random.random((100, 100))
    x2 = np.random.random((150, 175))

    # arbitrary extent - both images must have same extent if you want
    # them to be resampled into the same axes space
    extent = (0, 1, 0, 1)
    im1 = plt.imshow(x1, extent=extent)
    im2 = plt.imshow(x2, extent=extent)
    im2.set_visible(False)

    def toggle_images(event):
        'toggle the visible state of the two images'
        # if event.key != 't':
        #     return
        # print(ord(event.key))
        if event.key != ' ':  # could be any key; and why not use a block?
            return
        b1 = im1.get_visible()
        b2 = im2.get_visible()
        im1.set_visible(not b1)
        im2.set_visible(not b2)
        plt.draw()

    plt.connect('key_press_event', toggle_images)

    plt.show()


def mixed():
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm

    df_master = pd.read_csv('C:/24hrs/df_master.csv', sep=';')
    df_master_V = df_master
    model = sm.MixedLM.from_formula("Weight ~ Time", data, groups=data["Filter"])
    result = model.fit()
    print(result.summary())
    return result


def make_mixed_model_data():
    points = 80
    np.random.seed(1234)
    d = {'A': np.random.randn(points),
         'B': np.random.randn(points),
         'C': np.random.randn(points),
         'Ran': np.random.randint(0, 3, points),
         'Dep': 0}
    df = pd.DataFrame(d)
    df['Dep'] = 17 + 1*df.A + 2*df.B + 4*df.C + 5*(df.Ran-1) + 1*np.random.randn(len(df))
    categories = ['X', 'Y', 'Z']
    df['Ran'] = [categories[r] for r in df['Ran']]

    df_model = df[0:int(3*points/4)]
    df_test = df[len(df_model):]
    return df_model, df_test


def try_statsmodel(from_formula=True):
    """
    Scratch fn to test statsmodels pkg's mixed-model facilities (and my understanding of them).
    :return: 
    """
    df_model, df_test = make_mixed_model_data()  # df_test for out-of-sample testing

    if from_formula:
        # ----- Using statsmodels via formula (to get an intercept included in fit):
        fml = 'Dep ~ A + B + C'  # intercept should be implied (?) as in R::lme4
        model = sm.MixedLM.from_formula(fml, groups=df_model['Ran'], data=df_model)
    else:
        # ----- Using statsmodels MixedLM class directly, feeding it 3 pandas Dataframes:
        model = sm.MixedLM(endog=df_model['Dep'],           # endog  = Dependent variable
                           exog=df_model[['A', 'B', 'C']],  # exog   = Independent variables
                           groups=df_model['Ran'])          # groups = RE var(s) (here, as strings)
    """
    model = an instance of class MixedLM. It should have the following attributes
        (most are *NOT* methods, as claimed in the documentation--sigh).
        In the following: FE=fixed effects, RE=random effects. 
    ------
    model.data.endog = the input dep (Y) var values [ndarray]
              .exog = the input indep (X) FE-only var values [ndarray]
              .orig_endog = .endog as a Pandas Series
              .orig_exog = .exog as a Pandas Series
              .param_names = variable names, FE then RE [list of str]
         .group_labels = RE categories [list]
         .groups = RE values, as category values [Series]
         .row_indices = observation numbers in each RE group [dict(RE category: list of ints)]    
    """
    fit = model.fit()
    print(fit.summary())
    """
    fit = an instance of class MixedLMResults. It should have the following attributes
        (most are *NOT* methods, as claimed in the documentation--sigh).
        In the following: FE=fixed effects, RE=random effects. 
    ------
    fit.bse = std errors of FE & RE coefficients [Series, index = var name]
       .converged = True if fit converged [bool]
       .cov_re = covariance matrix of RE values [1x1 DataFrame]
       .fe_params = FE coefficient values from fit [Series, index = var name]
       .fittedvalues = fit values of dep. var, includes RE contribs [Series, length=# input pts]
       .llf = likelihood [float]
       .model = [object of MixedLM class]
       .nobs = number of observations included in fit [int]
       .params = FE then RE coefficient values from fit (RE is variance) [Series, index=var name]
       .params_object = [MixedLMParams object]
       .predict(exog=[DataFrame]) --> returns predicted Y  **ONLY INCLUDES FIXED-EFFECT CONTRIBS!**
       .random_effects = fit values of RE categories [dict(var name: float)]
       .resid = residuals of fit (Y_input - Y_fit) [Series, length = # input pts]
       .summary() -->  yields a summary table etc of fit [str?]
       .tvalues = t- (Z-) values of FE & RE coefficients from fit [Series, index = var name]   
    """

    # Get predicted values on fixed effects only (per statsmodels' weird def. of 'predicted'):
    fixed_effect_inputs = df_model[['A', 'B', 'C']]  # 1 col per fixed effect variable
    predicted_on_fixed_only = fit.predict(exog=fixed_effect_inputs)  # Series, length= # input pts

    # Get separate RE contib to pred values; for now allow only 1 random effect variable:
    random_effect_inputs = pd.DataFrame(df_model['Ran'])  # 1 col 'Ran' = the random-effect var
    random_effect_coeffs = pd.DataFrame(fit.random_effects).transpose()  # 1 row per group
    predicted_on_random_only = pd.merge(random_effect_inputs, random_effect_coeffs,
                                        left_on='Ran', right_index=True,
                                        how='left', sort=False)['groups']  # Series (a left-join)

    # Compare: constructed prediction values vs. expected:
    constructed_prediction = predicted_on_fixed_only + predicted_on_random_only
    expected_prediction = fit.fittedvalues
    pred_diff = constructed_prediction - expected_prediction
    print('max prediction diff = ', max(abs(pred_diff)))  # good, this is zero.

    # Compare: fitted values - input values (dep var)
    input_y = df_model['Dep']
    resid_fv = input_y - fit.fittedvalues
    resid_model = fit.resid
    df_resid = pd.DataFrame({'FromFV': resid_fv, 'FromModel': resid_model})
    max_diff = max(abs(resid_fv-resid_model))
    print('max resid diff = ', max_diff)  # good, they do match.

    # Out-of-sample prediction test:
    from math import sqrt
    rms_resid_in_sample = sqrt(sum(fit.resid*fit.resid)/(fit.nobs-1))
    pred_out_of_sample_fe = fit.predict(exog=df_test[['A', 'B', 'C']])
    re_inputs_out_of_sample = pd.DataFrame(df_test['Ran'])
    pred_out_of_sample_re = pd.merge(re_inputs_out_of_sample,
                                     random_effect_coeffs,  # this from in-sample code above.
                                     left_on='Ran', right_index=True,
                                     how='left', sort=False)['groups']  # Series (a left-join)
    pred_out_of_sample = pred_out_of_sample_fe + pred_out_of_sample_re
    rms_resid_out_of_sample = sqrt(sum((pred_out_of_sample - df_test['Dep'])**2)/(fit.nobs-1))
    print()
    print('rms(in sample resids) = ', rms_resid_in_sample)
    print('rms(out of sample resids) = ', rms_resid_out_of_sample)
    print()

    # Print relevant attributes of MixedLMResults object (named here 'fit'):
    print(30*'-')
    print('Attribute dump:')
    print(30 * '-')
    print('fit.bse = ', fit.bse)
    print('fit.converged = ', fit.converged)
    print('fit.cov_re = ', fit.cov_re)
    print('fit.fe_params = ', fit.fe_params)
    print('fit.fittedvalues.head() = ', fit.fittedvalues.head())
    print('fit.llf = ', fit.llf)
    print('fit.nobs = ', fit.nobs)
    print('fit.params = ', fit.params)
    print('fit.random_effects = ', fit.random_effects)
    print('fit.resid.head() = ', fit.resid.head())
    print('fit.tvalues = ', fit.tvalues)
    print(30*'-')


def try_logging():
    import logging
    import time
    logging.basicConfig(filename='c:/24hrs/logging.txt', level=logging.INFO,
                        format='%(asctime)s Z  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    logging.Formatter.converter = time.gmtime
    logging.info('----- precheck() J:/Astro/Images/C14/20170505')
    logging.info('188 FITS files found.')
    logging.info('12 FITS without plate solutions.')



