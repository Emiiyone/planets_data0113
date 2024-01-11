# m37.py
# マイコン宇宙講座
# 3-7 惑星の視位置計算プログラム
import math
import lib


M2PI = 2.0 * math.pi

lib.PL[3] = 'Sun'

x = [[0 for i in range(10)] for j in range(4)]
f = [0, 0, 0, 0]
q = [0, 0, 0, 0]

print()
## std = input('DATE AND TIME(JST) ? ')
## dy, dt = std.split(',')
## dy = float(dy)
## dt = float(dt)

std = input('DATE AND TIME(JST) ? [YYYY/MM/DD],[HH:MM:SS] ')
date_input, time_input = std.split(',')

## 日付の分割と結合
year, month, day = date_input.split('/')
dy = int(year + month + day)

## 時刻の分割と結合
hour, minute, second = time_input.split(':')
dt = int(hour + minute + second)

## ここで dy と dt を mjd 関数に渡す
jd, yy, mm, dd, hh, ms, ss = lib.mjd(dy, dt)



print()
print('  Date= %5d 年 %2d 月 %2d 日  Time= %2d 時 %2d 分 %2.1f 秒\n' % (yy, mm, dd, hh, ms, ss))
print('         MJD=%12.2f ET\n' % (jd))
print(' Planets        R.A.  (Date)   Decl.         Distance')
print(' ------------------------------------------------------')
print('                h   m            。  ,              AU')

t1 = jd - 33281.92334
t1 = t1 * (2.737909288e-5 + 1.260132857e-17 * t1)
t2 = t1 * t1

ty = yy + (mm - 1) / 12.0 + dd / 365.0
if yy < 0:
    ty = yy - (1 - ((mm - 1) / 12.0 + dd / 365.0))

for pn in range(1, 10):
    e, m, p, n, i, a, rd = lib.mean_elements(pn, t1, t2)
    f, q = lib.eph_const(p, n, i, lib.K)

    ec = e
    mo = M2PI * (m / (M2PI) - int(m / (M2PI)))

    # ss, cc, ff = lib.kepler(mo, ec)
    ss, cc, ff, _ = lib.kepler(mo, ec)


    b = a * math.sqrt(1.0 - e * e)
    for n in range(1, 4):
        f[n] = a * f[n]
        q[n] = b * q[n]

    # 惑星座標の計算
    for n in range(1, 4):
        x[n][pn] = ff * f[n] + ss * q[n]

for n in range(1, 4):
    x[n][3] = - (x[n][3])

for pn in range(1, 10):
    if pn != 3:
        for n in range(1, 4):
            x[n][pn] = x[n][pn] + x[n][3]


## ここからループで取得していくために追加
planets_data = []

# 惑星の赤経・赤緯の計算
for pn in range(1, 10):
    cc = x[1][pn]
    ss = x[2][pn]

    tt = lib.quadrant(ss, cc)

    ra = tt
    cc = x[1][pn] / math.cos(ra)
    ss = x[3][pn]

    dc = math.atan(ss / cc)

    r1 = ra
    d1 = dc

    ra1, dc1 = lib.precession(r1, d1, ty, 0, 0, lib.K)

    if ra1 > M2PI:
        ra1 = ra1 - M2PI
    if ra1 < 0:
        ra1 = ra1 + M2PI
    ra = ra1
    dc = dc1

    ra *= lib.K[3] / 15.0
    dc *= lib.K[3]
    ds = math.sqrt(x[1][pn]**2 + x[2][pn]**2 + x[3][pn]**2)

    rh = int(ra)
    rm = 60.0 * (ra - rh)
    dh = int(abs(dc))
    dm = 60.0 * (abs(dc) - dh)
    dh = lib.sgn(dc) * dh

    ## ここから追加
    degrees = rh * 15 + rm * 0.25 + (ra - rh - rm/60) * 15
    degree_part = int(degrees)
    minute_part = int((degrees - degree_part) * 60)
    ## ここまで
    ## print(' %-7s      %02d  %05.2f       %+03d  %04.1f       %8.5f' % (lib.PL[pn], rh, rm, dh, dm, ds))を変更
    print(' %-7s      %3d度%02d分' % (lib.PL[pn], degree_part, minute_part))

    print(' ------------------------------------------------------')
    print()

    ## jsonのためのコード追加

    import datetime
    import os
    import json

    ## JSONファイルに保存するための天体データをリスト化

    ## 必要な各天体のデータ
    planet_data = {
        "name": lib.PL[pn],
        "deg": degree_part,  # 各天体の計算結果（度数）を格納する
        "min": minute_part,  # 各天体の計算結果（分）を格納する
    }
    planets_data.append(planet_data)

    ## JSONファイルを保存するディレクトリの準備
    current_dir = os.path.dirname(os.path.abspath(__file__))
    save_dir = os.path.join(current_dir, "from_m37")
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    ## タイムスタンプ付きのファイル名を生成
    current_time = datetime.datetime.now()
    filename = current_time.strftime("planets_data_%Y%m%d_%H%M%S.json")

    ## ファイルの完全なパス
    full_path = os.path.join(save_dir, filename)

    ## データをファイルに保存
    with open(full_path, 'w') as file:
        json.dump(planets_data, file, indent=4)


