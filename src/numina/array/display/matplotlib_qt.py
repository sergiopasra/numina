# import matplotlib
# matplotlib.use('Qt5Agg')


def _pyplot():
    """Import matplotlib.pyplot, configured, the first time it is needed

    Importing pyplot takes some tenths of a second, and the modules that
    import this one only plot in some functions.
    """
    import matplotlib.pyplot as plt

    if "plt" not in globals():
        plt.rcParams.update({"figure.max_open_warning": 0})  # avoid warning
        globals()["plt"] = plt
    return plt


def __getattr__(name):
    # plt is imported on first access
    if name == "plt":
        return _pyplot()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def set_window_geometry(geometry):
    """Set window geometry.

    Parameters
    ==========
    geometry : tuple (4 integers) or None
        x, y, dx, dy values employed to set the Qt backend geometry.

    """

    if geometry is not None:
        x_geom, y_geom, dx_geom, dy_geom = geometry
        mngr = _pyplot().get_current_fig_manager()
        if "window" in dir(mngr):
            try:
                mngr.window.setGeometry(x_geom, y_geom, dx_geom, dy_geom)
            except AttributeError:
                pass
            else:
                pass
