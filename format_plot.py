def format_plot():
    import matplotlib as mpl

    # Use TeX for fonts
    # mpl.rcParams['text.usetex'] = True
    # Font settings for presentation slide graphics
    mpl.rcParams['font.serif'] = "Times"
    mpl.rcParams['font.family'] = "serif"
    mpl.rcParams['mathtext.fontset'] = "cm"
    mpl.rcParams['font.size'] = 18

    return
