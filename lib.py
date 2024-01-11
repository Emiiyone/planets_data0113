#  マイコン宇宙講座　サブルーチンライブラリ
# Python版
import math


# 月に対応させるためT[0]には0を入れてある。
T = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

# N-BASICの配列要素に対応
K = [0, 0.01720209895, 0, 180.0 / math.pi, 0.91743695, 0.39788118,
     (180.0 / math.pi) * 3600.0]

# 惑星名
PL = ['', 'MERCURY', 'VENUS', 'EARTH', 'MARS', 'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE', 'PLUTO']


def fno(x):
    return int(x + 0.5)


def fna(x):
    return int(x * 10.0 + 0.5) / 10.0


def fnb(x):
    return int(x * 100.0 + 0.5) / 100.0


def fnc(x):
    return int(x * 1000.0 + 0.5) / 1000.0


def fnr(x):
    return 2.0 * math.pi * (x / 360.0 - int(x / 360.0))


# N-BASICのSGN関数
def sgn(n):
    if n < 0:
        x = -1
    else:
        x = 1
    return x


# N-BASICのTAB関数
def tab(n):
    s = ''
    for i in range(1, n):
        s += ' '

    return s


# ユリウス日から年月日を求める
# 天文計算入門　63頁 (16.6)より
#
# 戻り値はyのみ、つまり年しか返していないが、mは月、dは日なので、
# 西暦年月日 = y年m月d日（UT)
# とすることができるのでreturn y, m, dとすれば年月日を返すこと
# ができる。
def get_dtime(jd, T):
    jd += 2400000.5

    a = int(jd + 68569.5)
    b = int(a / 36524.25)
    c = a - int(36524.25 * b + 0.75)
    e = int((c + 1) / 365.25025)
    f = c - int(365.25 * e) + 31
    g = int(f / 30.59)
    d = f - int(30.59 * g) + (jd + 0.5) - int(jd + 0.5)  # 日
    h = int(g / 11)
    m = g - 12 * h + 2  # 月
    y = 100 * (b - 49) + e + h  # 年

    # うるう年の判定
    if (y % 4 == 0 and y % 100 != 0) or y % 400 == 0:
        T[2] = 29

    # ユリウス日から年月日を求めるのでdの値がその月の日数を
    # 大きく超えることはなく、せいぜい＋1日程度だと思われる
    # 12月32日になったときの処理
    if m == 12:
        if d > T[m]:
            y += 1
            m = (m + 1) - 12
            d = 1
    else:
        # 12月以外の月における処理
        if d > T[m]:
            m += 1
            d = 1

    T[2] = 28

    return y, m, d


# MJD
def mjd(dy, dt):
    s1 = sgn(dy)
    dy = abs(dy)

    yy = int(dy / 10000.0)
    mm = int(dy / 100.0) % 100
    dd = dy % 100
    hh = int(dt / 10000.0)
    ms = int(dt / 100.0) % 100
    ss = dt % 100.0

    if s1 < 0:
        yy *= s1
    else:
        if yy < 100:
            yy += 1900

    jd = julian(yy, mm, dd)

    # ユリウス日に時刻を小数点にした値を加えてJSTからUTに変更
    jd += hh / 24.0 + ms / 1440.0 + ss / 86400.0
    jd -= 0.375

    # ユリウス日(JD)を修正ユリウス日(MJD)にする
    jd -= 2400000.5

    return jd, yy, mm, dd, hh, ms, ss


# JULIAN
def julian(yy, mm, dd):
    p5 = yy + (mm - 1.0) / 12.0 + dd / 365.25
    if mm < 3:
        mm += 12.0
        yy -= 1.0

    # 天文計算入門　長谷川一郎著　恒星社より引用
    # ３つの式をif文でネストすると警告がでるため、ネストしないで記述
    # ただし、記述する順番に気をつけないと正しい計算結果にならないので注意
    # また、マイコン宇宙講座の式と少し異なる。この式はユリウス日(JD)を求める式
    # ユリウス日(JD)と修正ユリウス日(MJD)の関係は以下のようになる
    #
    # MJD = JD - 2400000.5

    # yy < 0
    jd = sgn(yy) * int(abs(yy) * 365.25) + int(30.59 * (mm - 2)) + dd + 1721085.5

    # yy <= 1582.77
    if p5 >= 0:
        jd = int(yy * 365.25) + int(30.59 * (mm - 2)) + dd + 1721086.5
    # yy >= 1582.78
    # 天文計算においては、通常、この式を使用する
    if p5 >= 1582.78:
        jd = int(yy * 365.25) + int(yy / 400.0) - int(yy / 100.0) + int(30.59 * (mm - 2)) + dd + 1721088.5

    # N-BASICでは変数がすべてグローバル変数となるため、変更した場合はもとに戻す処理が必要
    # PHPでは引数を参照渡しする以外は不要
    yy *= sgn(yy)
    if mm > 12:
        mm -= 12.0
        yy += 1.0

    return jd


# JDATE
def jdate(jd, T):
    jj = jd
    yy = int(2.7379093e-3 * jd + 1858.877)
    mm = 1.0
    dd = 0.0

    jd = julian(yy, mm, dd) - 2400000.5
    r2 = jj - jd

    # うるう年の判定
    if (yy % 4 == 0 and yy % 100 != 0) or yy % 400 == 0:
        T[2] = 29

    r1 = 0
    m = 1
    while m < 13:
        if int(r2) - r1 - T[m] <= 0:
            break
        r1 += T[m]
        m += 1

    mm = m
    dd = r2 - r1
    T[2] = 28
    jd = jj

    if mm == 13:
        yy += 1.0
        mm -= 12.0

    return yy, mm, dd


# KOUSEI-JI
# 1950.0分点に準拠した地方恒星時
def siderealtime_1950(r1, r2):
    r3 = 0.6705199 + 1.002737803 * (r1 - 40000.0) + r2 / (2.0 * math.pi)
    r3 = r3 - int(r3)
    r3 = 2.0 * math.pi * r3

    return r3


# KOUSEI-JI
# 平均春分点に準拠した地方恒星時
def siderealtime_date(r1, r2):

    r3 = 0.671262 + 1.002737909 * (r1 - 40000.0) + r2 / (2.0 * math.pi)
    r3 = r3 - int(r3)
    r3 = 2.0 * math.pi * r3

    return r3


# SAISA
def precession(r1, d1, ty, p1, p2, K):
    p1 /= K[6]
    p2 /= K[6]
    ep = 1950.0
    t0 = (ep - 1900.0) / 100.0
    t1 = (ty - ep) / 100.0
    t2 = t1 * t1
    t3 = t2 * t1

    # 天体の位置計算 長沢工著 53頁　式3-11より
    x1 = (2304.25 + 1.396 * t0) * t1 + 0.302 * t2 + 0.018 * t3
    x2 = x1 + 0.791 * t2 + 0.001 * t3
    x3 = (2004.682 - 0.853 * t0) * t1 - 0.426 * t2 - 0.042 * t3

    x1 /= K[6]
    x2 /= K[6]
    x3 /= K[6]

    r2 = r1 + p1 * t1
    d2 = d1 + p2 * t1

    ss = math.cos(d2) * math.sin(r2 + x1)
    cc = math.cos(x3) * math.cos(d2) * math.cos(r2 + x1) - math.sin(x3) * math.sin(d2)

    r3 = quadrant(ss, cc)

    cc = ss / math.sin(r3)
    ss = math.cos(x3) * math.sin(d2) + math.sin(x3) * math.cos(d2) * math.cos(r2 + x1)
    dc = math.atan(ss / cc)
    ra = r3 + x2

    return ra, dc


# SYOGEN
def quadrant(ss, cc):
    tt = math.atan(ss / cc)
    if cc < 0:
        tt += math.pi
    else:
        if ss < 0:
            tt += (2.0 * math.pi)

    return tt


# 9惑星の軌道要素
#
def mean_elements(pn, t1, t2):
    # 水星
    if pn == 1:
        e = 0.20562441 + 2.042e-5 * t1 - 3e-8 * t2
        m = -0.729180963 + 2608.731797 * t1 + 8.242e-8 * t2
        p = 0.505094601 + 4.98049095e-3 * t1 + 1.33324e-6 * t2
        n = 0.833197496 - 2.1919881e-3 * t1 - 1.57564e-6 * t2
        i = 0.122238982 - 1.0510761e-4 * t1 + 1.454e-8 * t2
        a = 0.3870986011
        rd = 2439.0

    # 金星
    if pn == 2:
        e = 6.79676e-3 - 4.773e-5 * t1 + 9e-8 * t2
        m = -0.848503209 + 1021.306597 * t1 + 2.259232e-5 * t2
        p = 0.95357925 + 4.99862298e-3 * t1 - 2.077911e-5 * t2
        n = 1.330454065 - 4.85584535e-3 * t1 - 1.78896e-6 * t2
        i = 0.059237299 - 1.798659e-5 * t1 - 5.6723e-7 * t2
        a = 0.7233316286
        rd = 6052.0

    # 地球
    if pn == 3:
        e = 0.01673012 - 4.192e-5 * t1 - 1.3e-7 * t2
        m = -0.036149841 + 628.288592 * t1 - 2.86525e-6 * t2
        p = -1.262568078 + 9.7864e-3 * t1 + 2.55497e-6 * t2
        n = 3.044140013 - 4.21225519e-3 * t1 + 2.0847e-7 * t2
        i = 2.2713521e-4 * t1 - 2.618e-7 * t2
        a = 1.00000023
        rd = 695989.0

    # 火星
    if pn == 4:
        e = 0.09335426 + 9.056e-5 * t1 - 7e-8 * t2
        m = -3.326348457 + 334.0465218 * t1 + 3.08342e-6 * t2
        p = 4.991059739 + 0.01288077229 * t1 + 7.98973e-6 * t2
        n = 0.858200646 - 5.14920611e-3 * t1 - 1.107314e-5 * t2
        i = 0.032288058 - 1.4539562e-4 * t1 - 3.9755e-7 * t2
        a = 1.523688174
        rd = 3397.0

    # 木星
    if pn == 5:
        e = 0.04827062 + 4.7756e-5 * t1 + 2.2676e-5 * t2
        m = 5.28677165 + 52.9684478 * t1 + 9.996858e-5 * t2
        p = -1.50945736 - 1.7239975e-4 * t1 + 1.885925e-5 * t2
        n = 1.741507757 + 3.204618e-5 * t1 + 9.40539e-6 * t2
        i = 0.0228307 + 7.8055e-7 * t1 + 3.6991e-7 * t2
        a = 5.202833481
        rd = 71398.0

    # 土星
    if pn == 6:
        e = 0.05604508 - 2.5595e-5 * t1 - 1.6172e-5 * t2
        m = 1.165262433 + 21.328205912 * t1 - 5.4580324e-4 * t2
        p = -0.38321152 + 4.2707237e-4 * t1 + 2.6548394e-4 * t2
        n = 1.980742073 + 3.005845e-5 * t1 - 5.657776e-5 * t2
        i = 0.043422822 + 8.69756e-6 * t1 + 3.56623e-6 * t2
        a = 9.538762055
        rd = 60000.0

    # 天王星
    if pn == 7:
        e = 0.04613734 - 4.8118e-5 * t1 + 1.5396e-5 * t2
        m = -1.273611261 + 7.479637598 * t1 + 8.2743151e-4 * t2
        p = 1.716582467 - 2.37485982e-3 * t1 - 8.1424458e-4 * t2
        n = 1.28641873 + 6.429559e-4 * t1 + 3.97547e-6 * t2
        i = 0.013501673 - 1.72933e-5 * t1 - 8.7412e-7 * t2
        a = 19.19139128
        rd = 25400.0

    # 海王星
    if pn == 8:
        e = 9.71449e-3 + 1.095407e-3 * t1 + 3.62034e-4 * t2
        m = 2.724754505 + 3.994465294 * t1 + 0.048355704 * t2
        p = -1.62194737 - 0.18123685121 * t1 - 0.048352097 * t2
        n = 2.290559396 + 4.467073e-5 * t1 - 1.844231e-5 * t2
        i = 0.03096505 - 3.001e-6 * t1 + 3.6216e-7 * t2
        a = 30.06106906
        rd = 24300.0

    # 冥王星
    if pn == 9:
        e = 0.24824802 + 4.97082e-4 * t1 + 5.63208e-4 * t2
        m = -0.999328454 + 2.553487672 * t1 + 6.86515565e-3 * t2
        p = 1.977072713 - 0.01828324506 * t1 - 6.76126008e-3 * t2
        n = 1.913508742 + 1.027805e-5 * t1 + 5.6820159e-5 * t2
        i = 0.29929226 + 1.5198909e-04 * t1 + 5.269925e-5 * t2
        a = 39.52940243
        rd = 1500.0

    return e, m, p, n, i, a, rd


# KEPLER
def kepler(mo, ec):
    a1 = 0.0

    # Pythonではgoto文はサポートされていないので
    # while文で無限ループ処理してbreakでループを
    # 抜ける
    while True:
        a2 = ec * math.sin(a1 + mo)
        if abs(a2 - a1) > 1.0e-5:
            a1 = a2
        else:
            a1 = a2 + mo
            break

    ss = math.sin(a1)
    cc = math.cos(a1)
    ff = cc - ec

    return ss, cc, ff, a1


# EPH CONST
def eph_const(pe, nd, ic, K):
    f = [0, 0, 0, 0]
    q = [0, 0, 0, 0]

    r1 = math.sin(pe)
    r2 = math.sin(nd)
    r3 = math.sin(ic)
    r4 = math.cos(pe)
    r5 = math.cos(nd)
    r6 = math.cos(ic)

    f[1] = r5 * r4 - r1 * r6 * r2
    q[1] = -r1 * r5 - r4 * r6 * r2
    r7 = r6 * r5 * K[4]
    r8 = r2 * K[4]
    r9 = r3 * K[5]
    f[2] = r1 * r7 + r4 * r8 - r1 * r9
    q[2] = r4 * r7 - r1 * r8 - r4 * r9
    r7 = r3 * K[4]
    r8 = r6 * r5 * K[5]
    r9 = r2 * K[5]
    f[3] = r1 * r7 + r1 * r8 + r4 * r9
    q[3] = r4 * r7 + r4 * r8 - r1 * r9

    return f, q


# KOUDO ZAHYO
def ecliptic_coordinates(ra, dc, K):
    lx = math.cos(dc) * math.cos(ra)
    ly = math.sin(dc) * K[5] + math.cos(dc) * math.sin(ra) * K[4]
    lz = math.sin(dc) * K[4] - math.cos(dc) * math.sin(ra) * K[5]

    ss = ly
    cc = lx

    lp = quadrant(ss, cc)

    ss = lz
    cc = lx / math.cos(lp)
    bp = math.atan(ss / cc)

    return lp, bp


# KODO_HOUI
def altitude_direction(jd, lg, la, ra, dc):
    r1 = jd
    r2 = lg

    r3 = siderealtime_date(r1, r2)
    s1 = r3 - ra
    cc = -math.sin(dc) * math.cos(la) + math.cos(dc) * math.sin(la) * math.cos(s1)
    ss = math.cos(dc) * math.sin(s1)

    tt = quadrant(ss, cc)

    al = tt
    cc = cc / math.cos(al)
    ss = math.sin(dc) * math.sin(la) + math.cos(dc) * math.cos(la) * math.cos(s1)
    hi = math.atan(ss / cc)

    return al, hi


# 出没時刻および方位角の計算
def appear(jd, lg, la, ra, dc, ds, rd, sg, hd, K):
    if ds > 0:
        pp = math.atan(4.26e-5 / ds)
        rd = math.atan(rd / (1.5e+8 * ds))
    else:
        rd = 0.0
        pp = 0.0

    rr = 0.0102
    cc = -math.tan(dc) * math.tan(la)
    ss = math.sqrt(1 - cc * cc)
    h0 = math.atan(ss / cc)
    if h0 < 0:
        h0 = h0 + math.pi

    dh = (rr + rd - pp) / (math.cos(la) * math.cos(dc) * math.sin(h0))
    r1 = jd + sg
    r2 = lg

    r3 = siderealtime_date(r1, r2)

    s1 = r3

    if hd == 0:
        s2 = ra - h0 - dh
        ta = s2 - s1
        if ta < 0:
            ta = ta + 2.0 * math.pi
        ta = ta / 1.002738
        ta = K[3] * ta / 15.0
        if ta > 24:
            ta = ta - 24.0

        return ta, 0

    if hd == 1:
        s3 = ra + h0 + dh
        td = s3 - s1
        if td < 0:
            td = td + 2.0 * math.pi
        td = td / 1.002738
        td = K[3] * td / 15.0
        if td > 24:
            td = td - 24

        cc = math.sin(dc) / math.cos(la)
        ss = -math.cos(dc) * math.sin(h0)
        al = abs(math.atan(ss / cc))
        if dc < 0:
            al = -al + math.pi

        return td, al


# 地心緯度
def geocentric_latitude(r1, r2):
    return math.atan((0.9933055 + 1.1e-9 * r2) * math.tan(r1))


# 地心距離
def geocentric_distance(r1):
    return 6378.16 * (0.99832707 + 1.67644e-3 * math.cos(2.0 * r1) - 3.52e-6 * math.cos(4.0 * r1))


# 地理緯度
def geography_latitude(r1):
    r1 = math.atan(math.tan(r1) / 0.9933055)
    r2 = geocentric_distance(r1)

    return r2


# 軌道要素の表示
def print_elements(td, pe, nd, ic, ec, q, ep, T, K):
    jd = td
    yy, mm, dd = jdate(jd, T)
    a = q / (1.0 - ec)
    pd = 0.00995198 * pow(a, 1.5)
    no = 2.0 * math.pi * 86400.0 / pd
    mo = no * (td - ep)
    r1 = a / 6378.16
    ww = 0.174 * (2.0 - 2.5 * math.sin(ic)**2) / pow(r1, 3.5)
    nn = -0.174 * math.cos(ic) / pow(r1, 3.5)
    uu = no + ww
    r3 = 2.0 * math.pi * 86400.0 / uu
    if abs(mo * K[3]) > 0.1:
        print('    Epoch =  %4d  %2d  %6.3f ET   Mo = %9.3f' % (yy, mm, dd, mo))
    print('     T    =  %4d  %2d  %6.3f ET' % (yy, mm, dd))
    print('     Peri.=   %7.3f               e =     %7.5f' % (pe * K[3], ec))
    print('     Node =   %7.3f (Date)        a =%9.2f (Km)' % (nd * K[3], a))
    print('     Inc. =   %7.3f               n = %12.6f' % (ic * K[3], no * K[3]))
    print('        q =%9.2f (Km)           P = %8.2f 分' % (q, pd / 60.0))
    print('       ww =  %7.2f      Node Period = %8.2f 分' % (ww * K[3], r3 / 60.0))
    print('       nn =  %7.2f' % (nn * K[3]))
    print('\n')

    return a, pe, no, pd, ww, nn


# KOUDO ZAHYO from XYZ
def ecliptic_coordinates_xyz(x, p, K):
    xe = x[1][p]
    ye = x[2][p] * K[4] + x[3][p] * K[5]
    ze = -x[2][p] * K[5] + x[3][p] * K[4]

    ss = ye
    cc = xe
    tt = quadrant(ss, cc)
    lp = tt
    cc = ye / math.sin(lp)
    ss = ze
    bp = math.atan(ss / cc)
    r0 = math.sqrt(xe * xe + ye * ye + ze * ze)

    return lp, bp, r0


def print_elements_display(td, pe, ec, nd, ic, no, q, K, T):
    jd = td
    yy, mm, dd = jdate(jd, T)
    a = q / (1.0 - ec)
    pd = pow(a, 1.5)
    print('     T    =  %4d  %2d  %6.3f ET' % (yy, mm, dd))
    print('     Peri.= %7.3f               e = %9.5f' % (pe * K[3], ec))
    print('     Node = %7.3f (1950)        a = %9.5f' % (nd * K[3], a))
    print('     Inc. = %7.3f               n =  %9.6f' % (ic * K[3], no * K[3]))
    print('        q =  %8.5f (AU)        P =   %4.2f 年' % (q, pd))

    return a


def argument(ta):
    tb = ta * ta
    a = 296.104608 + 477000.0 * ta + 198.849108 * ta + 9.192e-3 * tb
    a = fnr(a)
    b = 11.250889 + 483120.0 * ta + 82.02515 * ta - 3.211e-3 * tb
    b = fnr(b)
    c = 270.434164 + 480960.0 * ta + 307.883142 * ta - 1.133e-3 * tb
    c = fnr(c)
    d = 350.737486 + 444960.0 * ta + 307.114217 * ta - 1.436e-3 * tb
    d = fnr(d)
    e = 98.998753 + 35640.0 * ta + 359.372886 * ta
    e = fnr(e)
    g = 358.475833 + 35999.04975 * ta - 1.5e-4 * tb
    g = fnr(g)
    j = 225.444651 + 2880.0 * ta + 154.906654 * ta
    j = fnr(j)
    l = 279.696678 + 36000.76892 * ta + 3.03e-4 * tb
    l = fnr(l)
    m = 319.529425 + 19080.0 * ta + 59.8585 * ta + 1.81e-4 * tb
    m = fnr(m)
    n = 259.183275 - 1800.0 * ta - 134.142008 * ta + 2.078e-3 * tb
    n = fnr(n)
    v = 212.603219 + 58320.0 * ta + 197.803875 * ta + 1.286e-3 * tb
    v = fnr(v)
    w = 342.767053 + 58320.0 * ta + 199.211911 * ta + 3.1e-4 * tb
    w = fnr(w)

    return a, b, c, d, e, g, j, l, m, n, v, w


def moon(ta, a, b, d, e, n, g, w):
    # MOON-THETA
    mx = -1.01e-3 * math.sin(a + b - 4 * d)
    mx = mx - 1.02e-3 * math.sin(a - b - 4 * d - n)
    mx = mx - 1.03e-3 * ta * math.sin(a - b - n)
    mx = mx - 1.07e-3 * math.sin(a - g - b - 2 * d - n)
    mx = mx - 1.21e-3 * math.sin(2.0 * a - b - 4 * d - n)
    mx = mx + 1.30e-3 * math.sin(3 * a + b + n)
    mx = mx - 1.31e-3 * math.sin(a + b - n)
    mx = mx + 1.36e-3 * math.sin(a + b - d + n)
    mx = mx - 1.45e-3 * math.sin(g + b)
    mx = mx - 1.49e-3 * math.sin(a + g - b - 2 * d)
    mx = mx + 1.57e-3 * math.sin(g - b + d - n)
    mx = mx - 1.59e-3 * math.sin(g - b)
    mx = mx + 1.84e-3 * math.sin(a - g + b - 2 * d + n)
    mx = mx - 1.94e-3 * math.sin(b - 2 * d - n)
    mx = mx - 1.96e-3 * math.sin(g - b + 2 * d - n)
    mx = mx + 2.00e-3 * math.sin(b - d)
    mx = mx - 2.05e-3 * math.sin(a + g - b)
    mx = mx + 2.35e-3 * math.sin(a - g - b)
    mx = mx + 2.46e-3 * math.sin(a - 3 * b - n)
    mx = mx - 2.62e-3 * math.sin(2 * a + b - 2 * d)
    mx = mx - 2.83e-3 * math.sin(a + g + b - 2 * d)
    mx = mx - 3.39e-3 * math.sin(g - b - 2 * d - n)
    mx = mx + 3.45e-3 * math.sin(a - b + n)
    mx = mx - 3.47e-3 * math.sin(g - b + 2 * d)
    mx = mx - 3.83e-3 * math.sin(b + d + n)
    mx = mx - 4.11e-3 * math.sin(a + g + b + n)
    mx = mx - 4.42e-3 * math.sin(2 * a - b - 2 * d - n)
    mx = mx + 4.49e-3 * math.sin(a - b + 2 * d)
    mx = mx - 4.56e-3 * math.sin(3 * b - 2 * d + n)
    mx = mx + 4.66e-3 * math.sin(a + b + 2 * d + n)
    mx = mx + 4.9e-3 * math.sin(2 * a - b)
    mx = mx + 5.61e-3 * math.sin(2 * a + b)
    mx = mx + 5.64e-3 * math.sin(a - g + b + n)
    mx = mx - 6.38e-3 * math.sin(a + g - b - n)
    mx = mx - 7.13e-3 * math.sin(a + g - b - 2 * d - n)
    mx = mx - 9.29e-3 * math.sin(g + b - 2 * d)
    mx = mx - 9.47e-3 * math.sin(2 * a - b - n)
    mx = mx + 9.65e-3 * math.sin(a - g - b - n)
    mx = mx + 9.70e-3 * math.sin(b + 2 * d)
    mx = mx + 0.01064 * math.sin(b - d + n)
    mx = mx - 0.0125 * ta * math.sin(b + n)
    mx = mx - 0.01434 * math.sin(g + b - 2 * d + n)
    mx = mx - 0.01652 * math.sin(a + g + b - 2 * d + n)
    mx = mx - 0.01868 * math.sin(2 * a + b - 2 * d + n)
    mx = mx + 0.02700 * math.sin(2 * a + b + n)
    mx = mx - 0.02994 * math.sin(a - b - 2 * d)
    mx = mx - 0.03759 * math.sin(g + b + n)
    mx = mx - 0.03982 * math.sin(g - b - n)
    mx = mx + 0.04732 * math.sin(b + 2 * d + n)
    mx = mx - 0.04771 * math.sin(b - n)
    mx = mx - 0.06505 * math.sin(a + b - 2 * d)
    mx = mx + 0.13622 * math.sin(a + b)
    mx = mx - 0.14511 * math.sin(a - b - 2 * d - n)
    mx = mx - 0.18354 * math.sin(b - 2 * d)
    mx = mx - 0.20017 * math.sin(b - 2 * d + n)
    mx = mx - 0.38899 * math.sin(a + b - 2 * d + n)
    mx = mx + 0.40248 * math.sin(a - b)
    mx = mx + 0.65973 * math.sin(a + b + n)
    mx = mx + 1.96763 * math.sin(a - b - n)
    mx = mx + 4.95372 * math.sin(b)
    mx = mx + 23.89684 * math.sin(b + n)

    # MOON-RHO
    my = 0.05491 * math.cos(2 * a + g)
    my = my + 0.06290 * math.cos(a + d)
    my = my - 0.06444 * math.cos(4 * d)
    my = my - 0.06652 * math.cos(2 * a - g)
    my = my - 0.07369 * math.cos(g - 4 * d)
    my = my + 0.08119 * math.cos(a - 3 * d)
    my = my - 0.09261 * math.cos(a + 4 * g)
    my = my + 0.10177 * math.cos(a - 2 * b + 2 * d)
    my = my + 0.10225 * math.cos(a + g + 2 * d)
    my = my - 0.10243 * math.cos(a + 2 * g - 2 * d)
    my = my - 0.12291 * math.cos(2 * b)
    my = my - 0.12291 * math.cos(2 * a - 2 * b)
    my = my - 0.12428 * math.cos(a + g - 4 * d)
    my = my - 0.14986 * math.cos(3 * a)
    my = my - 0.16070 * math.cos(a - g + 2 * d)
    my = my - 0.16949 * math.cos(a - d)
    my = my + 0.17697 * math.cos(a + 2 * b - 2 * d)
    my = my - 0.18815 * math.cos(2 * a - 4 * d)
    my = my - 0.19482 * math.cos(2 * g - 2 * d)
    my = my + 0.22383 * math.cos(2 * b - 2 * d)
    my = my + 0.22594 * math.cos(3 * a - 2 * d)
    my = my + 0.24454 * math.cos(2 * a + g - 2 * d)
    my = my - 0.31717 * math.cos(g + d)
    my = my - 0.36333 * math.cos(a - 4 * d)
    my = my + 0.47999 * math.cos(a - g - 2 * d)
    my = my + 0.63844 * math.cos(g + 2 * d)
    my = my + 0.86170 * math.cos(g)
    my = my + 1.50534 * math.cos(a - 2 * b)
    my = my - 1.67417 * math.cos(a + 2 * d)
    my = my + 1.99463 * math.cos(a + g)
    my = my + 2.07579 * math.cos(d)
    my = my - 2.45500 * math.cos(a - g)
    my = my - 2.74067 * math.cos(a + g - 2 * d)
    my = my - 3.83002 * math.cos(g - 2 * d)
    my = my - 5.37817 * math.cos(2 * a)
    my = my + 6.60763 * math.cos(2 * a - 2 * d)
    my = my - 53.97626 * math.cos(2 * d)
    my = my - 68.62152 * math.cos(a - 2 * d)
    my = my - 395.13669 * math.cos(a)
    my = my + 3649.33705

    # MOON-PHI
    mz = -1.00e-3 * math.sin(a - g - 2 * b - 2 * n)
    mz = mz - 1.00e-3 * math.sin(a + g - 4 * d)
    mz = mz + 1.00e-3 * math.sin(2 * a - g)
    mz = mz + 1.02e-3 * math.sin(a - g + 2 * d)
    mz = mz - 1.06e-3 * math.sin(2 * a - 2 * b - n)
    mz = mz - 1.06e-3 * math.sin(2 * a + n)
    mz = mz - 1.09e-3 * math.sin(a + 2 * b - 2 * d)
    mz = mz - 1.10e-3 * math.sin(2 * b - d + 2 * n)
    mz = mz + 1.12e-3 * math.sin(4 * d)
    mz = mz - 1.22e-3 * math.sin(2 * a - n)
    mz = mz - 1.22e-3 * math.sin(2 * a + 2 * b + n)
    mz = mz + 1.49e-3 * math.sin(g + 2 * b - 2 * d + 2 * n)
    mz = mz - 1.57e-3 * math.sin(2 * a - 4 * d)
    mz = mz + 1.71e-3 * math.sin(a + g + 2 * b - 2 * d + 2 * n)
    mz = mz - 1.75e-3 * math.sin(2 * a + g - 2 * d)
    mz = mz - 1.90e-3 * math.sin(2 * g - 2 * d)
    mz = mz + 1.93e-3 * math.cos(a + 16 * e - 18 * w)
    mz = mz + 1.94e-3 * math.sin(2 * a + 2 * b - 2 * d + 2 * n)
    mz = mz + 2.01e-3 * math.sin(g - 2 * d - n)
    mz = mz + 2.01e-3 * math.sin(g + 2 * b - 2 * d + n)
    mz = mz - 2.07e-3 * math.sin(a + 2 * g - 2 * d)
    mz = mz - 2.10e-3 * math.sin(2 * g)
    mz = mz - 2.13e-3 * math.sin(2 * d - n)
    mz = mz - 2.13e-3 * math.sin(2 * b + 2 * d + n)
    mz = mz - 2.15e-3 * math.sin(3 * a - 2 * d)
    mz = mz - 2.47e-3 * math.sin(a - 4 * d)
    mz = mz - 2.53e-3 * math.sin(a - 2 * b + 2 * d)
    mz = mz + 2.79e-3 * ta * math.sin(2 * b + 2 * n)
    mz = mz - 2.80e-3 * math.sin(2 * a + 2 * b + 2 * n)
    mz = mz + 3.12e-3 * math.sin(3 * a)
    mz = mz - 3.17e-3 * math.sin(a + 2 * b)
    mz = mz - 3.50e-3 * math.sin(a + 16 * e - 18 * w)
    mz = mz + 3.90e-3 * math.sin(g + 2 * b + 2 * n)
    mz = mz + 4.13e-3 * math.sin(g - 2 * b - 2 * n)
    mz = mz - 4.90e-3 * math.sin(2 * n)
    mz = mz - 4.91e-3 * math.sin(2 * b + 2 * d + 2 * n)
    mz = mz + 5.04e-3 * math.sin(g + d)
    mz = mz + 5.16e-3 * math.sin(a - d)
    mz = mz - 6.21e-3 * math.sin(g + 2 * d)
    mz = mz + 6.48e-3 * math.sin(a - 2 * b - 2 * d - n)
    mz = mz + 6.48e-3 * math.sin(a - 2 * d + n)
    mz = mz + 7.00e-3 * math.sin(a - g - 2 * d)
    mz = mz + 0.01122 * math.sin(a + 2 * d)
    mz = mz + 0.01410 * math.sin(a - 2 * d - n)
    mz = mz + 0.01410 * math.sin(a + 2 * b - 2 * d + n)
    mz = mz + 0.01424 * math.sin(a - 2 * b)
    mz = mz + 0.01506 * math.sin(a - 2 * b - 2 * d - 2 * n)
    mz = mz - 0.01567 * math.sin(2 * b - 2 * d)
    mz = mz + 0.02077 * math.sin(2 * b - 2 * d + 2 * n)
    mz = mz - 0.02527 * math.sin(a + g)
    mz = mz - 0.02952 * math.sin(a - n)
    mz = mz - 0.02952 * math.sin(a + 2 * b + n)
    mz = mz - 0.03487 * math.sin(d)
    mz = mz + 0.03684 * math.sin(a - g)
    mz = mz - 0.03983 * math.sin(2 * d + n)
    mz = mz + 0.03983 * math.sin(2 * b - 2 * d + n)
    mz = mz + 0.04037 * math.sin(a + 2 * b - 2 * d + 2 * n)
    mz = mz + 0.04221 * math.sin(2 * a)
    mz = mz - 0.04273 * math.sin(g - 2 * d)
    mz = mz - 0.05566 * math.sin(2 * a - 2 * d)
    mz = mz - 0.05697 * math.sin(a + g - 2 * d)
    mz = mz - 0.06846 * math.sin(a + 2 * b + 2 * n)
    mz = mz - 0.08724 * math.sin(a - 2 * b - n)
    mz = mz - 0.08724 * math.sin(a + n)
    mz = mz - 0.11463 * math.sin(2 * b)
    mz = mz - 0.18647 * math.sin(g)
    mz = mz - 0.20417 * math.sin(a - 2 * b - 2 * n)
    mz = mz + 0.59616 * math.sin(2 * d)
    mz = mz + 1.07142 * math.sin(n)
    mz = mz - 1.07447 * math.sin(2 * b + n)
    mz = mz - 1.28658 * math.sin(a - 2 * d)
    mz = mz - 2.47970 * math.sin(2 * b + 2 * n)
    mz = mz + 6.32962 * math.sin(a)

    return mx, my, mz


# POSITIONS
def positions(r0, r1, r2, r3):
    ss = r3 / math.sqrt(r2 - r1 * r1)
    cc = math.sqrt(1.0 - ss * ss)
    ra = r0 + math.atan(ss / cc)

    ss = r1 / math.sqrt(r2)
    cc = math.sqrt(1.0 - ss * ss)
    dc = math.atan(ss / cc)

    dl = math.sqrt(r2)

    return ra, dc, dl


# SUN
def sun(ta, c, g, j, l, m, n, v):
    # SUN-THETA
    sx = - 1.000e-5 * ta * math.sin(g + l)
    sx = sx - 1.000e-5 * math.cos(g - l - j)
    sx = sx - 1.400e-5 * math.sin(2 * g - l)
    sx = sx - 3.000e-5 * ta * math.sin(g - l)
    sx = sx - 3.900e-5 * math.sin(n - l)
    sx = sx - 4.000e-5 * math.cos(l)
    sx = sx + 4.200e-5 * math.sin(2 * g + l)
    sx = sx - 2.080e-4 * ta * math.sin(l)
    sx = sx + 3.334e-3 * math.sin(g + l)
    sx = sx + 9.999e-3 * math.sin(g - l)
    sx = sx + 0.39793 * math.sin(l)

    # SUN-RHO
    sy = 2.7e-5 * math.sin(2 * g - 2 * v)
    sy = sy - 3.3e-5 * math.sin(g - j)
    sy = sy + 8.4e-5 * ta * math.cos(g)
    sy = sy - 1.4e-4 * math.cos(2 * g)
    sy = sy - 0.033503 * math.cos(g)
    sy = sy + 1.000421

    # SUN-PHI
    sz = -1.700e-5 * math.cos(2 * g - 2 * v)
    sz = sz - 1.900e-5 * math.sin(g - v)
    sz = sz + 2.400e-5 * math.sin(4 * g - 8 * m + 3 * j)
    sz = sz - 2.500e-5 * math.cos(g - j)
    sz = sz + 3.000e-5 * math.sin(c - l)
    sz = sz + 4.600e-5 * ta * math.sin(2 * l)
    sz = sz + 6.801e-5 * math.sin(2 * g)
    sz = sz - 7.900e-5 * math.sin(n)
    sz = sz - 8.000e-5 * ta * math.sin(g)
    sz = sz - 9.500e-5
    sz = sz - 3.460e-4 * math.sin(g + 2 * l)
    sz = sz - 1.038e-3 * math.sin(g - 2 * l)
    sz = sz + 0.032116 * math.sin(g)
    sz = sz - 0.041295 * math.sin(2 * l)

    return sx, sy, sz


# NEARLY PARABOLIC
def nearly_parabolic(tp, t, q, ec, K):
    r1 = 1.0 + 9.0 * ec
    aa = math.sqrt(0.1 * r1)
    bb = 5.0 * (1.0 - ec) / r1
    cc = math.sqrt(5.0 * (1.0 + ec) / r1)

    b = 1.0
    a0 = 0.0
    while True:
        u = b * aa * K[1] * (tp - t) / (math.sqrt(2.0) * math.pow(q, 1.5))
        r2 = 1.0
        while True:
            v = (u + 2.0 * math.pow(r2, 3) / 3.0) / (1.0 + r2 * r2)
            if abs(v - r2) > 1.0e-6:
                r2 = v
            else:
                break
        a = bb * v * v
        a2 = a * a
        a3 = a2 * a
        b = 1.0 - 0.017142857 * a2 - 0.003809524 * a3
        if abs(a - a0) > 1.0e-6:
            a0 = a
        else:
            break

    c = 1.0 + 0.4 * a + 0.21714286 * a2 + 0.12495238 * a3
    d = 1.0 - a + 0.2 * a2 + 0.00571429 * a3
    qd = q * d
    v = cc * c * v

    cc = qd * (1.0 - math.pow(v, 2))
    ss = 2.0 * qd * v
    v = 2.0 * math.atan(v)

    return ss, cc, v
