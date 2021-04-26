import input_output as io


def sod_solve(zz, yy, nmin, nmax, h=None):
    r12 = 1.0 / 12.0

    # Decide whether integrating from:
    #    left to right ---> isgn = + 1
    # or right to left ---> isgn = - 1
    isgn = (nmax - nmin) // abs(nmax - nmin)

    # Run some test just to be conservative.
    if isgn == +1:
        if nmin <= 2 or nmax > len(zz):
            io.abort('sod_solve: nmin, nmax out of range')

    elif isgn == -1:
        if nmin >= (len(zz) - 1) or nmax < 1:
            io.abort('sod_solve: nmin, nmax out of range')

    else:
        io.abort('sod_solve: nmin, nmax out of range')

    # Initialise current second derivative of yy (d2y), ycur and yold.
    if h is not None:
        i = nmin - isgn - 1
        d2y = zz[i] * yy[i] + h[i]
        ycur = yy[i] - r12 * (zz[i] * yy[i] + h[i])
        i -= isgn
        yold = yy[i] - r12 * (zz[i] * yy[i] + h[i])

    else:
        i = nmin - isgn - 1
        d2y = zz[i] * yy[i]
        ycur = (1.0 - r12 * zz[i]) * yy[i]
        i -= isgn
        yold = (1.0 - r12 * zz[i] ) * yy[i]

    # Begin the integration loop
    if h is not None:
        for i in range(nmin-isgn, nmax+isgn, isgn):
            ynew  = 2.0 * ycur - yold + d2y
            yy[i] = (ynew + r12 * h[i]) / (1.0 - r12 * zz[i])
            d2y   = zz[i] * yy[i] + h[i]
            yold  = ycur
            ycur  = ynew

    else:
        for i in range(nmin-isgn, nmax+isgn, isgn):
            ynew  = 2.0 * ycur - yold + d2y
            yy[i] = ynew / (1.0 - r12 * zz[i])
            d2y   = zz[i] * yy[i]
            yold  = ycur
            ycur  = ynew

