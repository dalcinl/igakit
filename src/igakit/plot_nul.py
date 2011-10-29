class colors(object):
    black   = (0.00, 0.00, 0.00)
    grey    = (0.75, 0.75, 0.75)
    white   = (1.00, 1.00, 1.00)
    #
    red     = (1.00, 0.00, 0.00)
    green   = (0.00, 1.00, 0.00)
    blue    = (0.00, 0.00, 1.00)
    #
    yellow  = (1.00, 1.00, 0.00)
    cyan    = (0.00, 1.00, 1.00)
    magenta = (1.00, 0.00, 1.00)
    #
    orange  = (1.00, 0.50, 0.00)

class _namespace(object):
    def __getattr__(self, name):
        def func(*args, **kwargs): return None
        func.__name__ = name
        func.__doc__  = ''
        return func
_namespace = _namespace()

show = _namespace.show

figure = _namespace.figure
close = _namespace.close
save = _namespace.savefig

title = _namespace.title
xlabel = _namespace.xlabel
ylabel = _namespace.ylabel
zlabel = _namespace.zlabel
colorbar = _namespace.colorbar

points3d = _namespace.points3d
quiver3d = _namespace.quiver3d
line3d = _namespace.quiver3d
surf3d = _namespace.quiver3d

_resolution = { 1:64, 2:32 }
