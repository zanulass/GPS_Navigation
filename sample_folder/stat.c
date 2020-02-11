/*------------------------------------------------------------
 * STAT.c - Position Statistics.
 *------------------------------------------------------------*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*------------------------------------------------------------
 * 定数・構造体の定義
 *------------------------------------------------------------*/
/* 論理型 */
typedef	int	bool;
#define	TRUE	1
#define	FALSE	0

/* WGS-84定数 */
#define	PI		3.1415926535898		/* 円周率（IS-GPS-200） */
#define	C		2.99792458e8		/* 光速[m/s] */
#define	MUe		3.986005e14			/* 地球重力定数[m^3/s^2] */
#define	dOMEGAe	7.2921151467e-5		/* 地球自転角速度[rad/s] */
#define	Re		6378137.0			/* 地球半径[m] */
#define	Fe		(1.0/298.257223563)	/* 地球の扁平率 */

/* 角度の変換 */
#define	rad_to_deg(rad)		((rad)/PI*180.0)
#define	deg_to_rad(deg)		((deg)/180.0*PI)
#define	rad_to_sc(rad)		((rad)/PI)
#define	sc_to_rad(sc)		((sc)*PI)

/* 取り扱える行列の大きさ */
#define	MAX_N				16		/* 観測衛星数の上限 */
#define	MAX_M				4		/* 未知数の最大数 */
#define	MAX_PRN				32		/* 衛星番号の上限 */

/* 時間 */
#define	SECONDS_DAY			(3600L*24L)
#define	SECONDS_WEEK		(3600L*24L*7L)

/* 時刻を表す構造体 */
typedef struct {
	int		week;		/* 週番号 */
	double	sec;		/* 週初めからの経過時間[s] */
}	wtime;

/* 直交座標を表す構造体 */
typedef struct {
	double	x;			/* X座標[m] */
	double	y;			/* Y座標[m] */
	double	z;			/* Z座標[m] */
}	posxyz;
#define	SQ(x)			((x)*(x))
#define	DIST(a,b)		sqrt(SQ(a.x-b.x)+SQ(a.y-b.y)+SQ(a.z-b.z))

/* 経緯度を表す構造体 */
typedef struct {
	double	lat;		/* 緯度[rad] */
	double	lon;		/* 経度[rad] */
	double	hgt;		/* 高度（楕円体高）[m] */
}	posblh;

/* ENU座標を表す構造体 */
typedef struct {
	double	e;			/* East成分[m] */
	double	n;			/* North成分[m] */
	double	u;			/* Up成分[m] */
}	posenu;

/*------------------------------------------------------------
 * ファイルからの読込み
 *------------------------------------------------------------*/
/* 行バッファ */
#define	LINEBUF_LEN			256
static char	linebuf[LINEBUF_LEN],fieldbuf[LINEBUF_LEN];
static int	linepos			=0;

/* ファイルから行バッファに読み込む */
static bool read_line(FILE *fp)
{
	/* 行バッファをクリア */
	linepos		=0;
	linebuf[0]	='\0';

	/* ファイルから読み込む */
	return (fgets(linebuf,LINEBUF_LEN,fp)!=NULL);
}

/* 行バッファから文字列を取り出す（文字数指定、またはCSV形式） */
static char *get_field(int width)
{
	int		i;
	bool	quotef=TRUE;

	/* width==0ならCSV形式 */
	if (width<1) {
		width	=LINEBUF_LEN-1;
		quotef	=FALSE;
	}

	/* 指定された文字数あるいはコンマまでを読み取る */
	for(i=0;i<width && linebuf[linepos+i]!='\0';i++) {
		if (linebuf[linepos+i]=='\n') break;
		if (!quotef && linebuf[linepos+i]==',') {
			linepos++;
			break;
		}
		if (linebuf[linepos+i]=='\"') {
			quotef=!quotef;
			continue;
		}
		fieldbuf[i]=linebuf[linepos+i];
	}
	linepos		+=i;		/* 次の位置 */
	fieldbuf[i]	='\0';

	return fieldbuf;
}

/*------------------------------------------------------------
 * 座標変換
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * xyz_to_blh() - 直交座標から経緯度への変換
 *------------------------------------------------------------
 *  posblh xyz_to_blh(pos); 経緯度
 *    posxyz pos; 直交座標値
 *------------------------------------------------------------*/
posblh xyz_to_blh(posxyz pos)
{
	double	a,b,e,f,n,h,p,t,sint,cost;
	posblh	blh={0.0,0.0,-Re};

	/* 原点の場合 */
	if (pos.x==0.0 && pos.y==0.0 && pos.z==0.0) return blh;

	/* 楕円体のパラメータ */
	f	=Fe;				/* 扁平率 */
	a	=Re;				/* 長半径 */
	b	=a*(1.0-f);			/* 短半径 */
	e	=sqrt(f*(2.0-f));	/* 離心率 */

	/* 座標変換のためのパラメータ */
	h	=a*a-b*b;
	p	=sqrt(pos.x*pos.x+pos.y*pos.y);
	t	=atan2(pos.z*a,p*b);
	sint=sin(t);
	cost=cos(t);

	/* 経緯度への変換 */
	blh.lat =atan2(pos.z+h/b*sint*sint*sint,p-h/a*cost*cost*cost);
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	blh.lon	=atan2(pos.y,pos.x);
	blh.hgt	=(p/cos(blh.lat))-n;

	return blh;
}

/*------------------------------------------------------------
 * blh_to_xyz() - 経緯度から直交座標への変換
 *------------------------------------------------------------
 *  posxyz blh_to_xyz(blh); 直交座標値
 *    posblh blh; 経緯度
 *------------------------------------------------------------*/
posxyz blh_to_xyz(posblh blh)
{
	double	a,b,e,f,n;
	posxyz	pos;

	/* 楕円体のパラメータ */
	f	=Fe;				/* 扁平率 */
	a	=Re;				/* 長半径 */
	b	=a*(1.0-f);			/* 短半径 */
	e	=sqrt(f*(2.0-f));	/* 離心率 */

	/* 直交座標系への変換 */
	n		=a/sqrt(1.0-e*e*sin(blh.lat)*sin(blh.lat));
	pos.x	=(n+blh.hgt)*cos(blh.lat)*cos(blh.lon);
	pos.y	=(n+blh.hgt)*cos(blh.lat)*sin(blh.lon);
	pos.z	=(n*(1.0-e*e)+blh.hgt)*sin(blh.lat);

	return pos;
}

/*------------------------------------------------------------
 * xyz_to_enu() - 直交座標からENU座標への変換
 *------------------------------------------------------------
 *  posenu xyz_to_enu(pos,base); ENU座標値
 *    posxyz pos;  直交座標値
 *    posxyz base; 基準位置
 *------------------------------------------------------------*/
posenu xyz_to_enu(posxyz pos,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posenu	enu;

	/* 基準位置からの相対位置 */
	pos.x	-=base.x;
	pos.y	-=base.y;
	pos.z	-=base.z;

	/* 基準位置の経緯度 */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* 相対位置を回転させてENU座標に変換する */
	enu.e	=-pos.x*s1+pos.y*c1;
	enu.n	=-pos.x*c1*s2-pos.y*s1*s2+pos.z*c2;
	enu.u	=pos.x*c1*c2+pos.y*s1*c2+pos.z*s2;

	return enu;
}

/*------------------------------------------------------------
 * enu_to_xyz() - ENU座標から直交座標への変換
 *------------------------------------------------------------
 *  posxyz enu_to_xyz(enu,base); 直交座標値
 *    posenu enu;  ENU座標値
 *    posxyz base; 基準位置
 *------------------------------------------------------------*/
posxyz enu_to_xyz(posenu enu,posxyz base)
{
	double	s1,c1,s2,c2;
	posblh	blh;
	posxyz	pos;

	/* 基準位置の経緯度 */
	blh	=xyz_to_blh(base);
	s1	=sin(blh.lon);
	c1	=cos(blh.lon);
	s2	=sin(blh.lat);
	c2	=cos(blh.lat);

	/* ENU座標を回転させて相対位置に変換する */
	pos.x	=-enu.e*s1-enu.n*c1*s2+enu.u*c1*c2;
	pos.y	=enu.e*c1-enu.n*s1*s2+enu.u*s1*c2;
	pos.z	=enu.n*c2+enu.u*s2;

	/* 基準位置に加える */
	pos.x	+=base.x;
	pos.y	+=base.y;
	pos.z	+=base.z;

	return pos;
}

/*------------------------------------------------------------
 * elevation() - 仰角を求める
 *------------------------------------------------------------
 *  double elevation(sat,usr); 仰角[rad]
 *    posxyz sat; 衛星の位置
 *    posxyz usr; ユーザ位置
 *------------------------------------------------------------*/
double elevation(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU座標に変換して仰角を求める */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.u,sqrt(enu.e*enu.e+enu.n*enu.n));
}

/*------------------------------------------------------------
 * azimuth() - 方位角を求める
 *------------------------------------------------------------
 *  double azimuth(sat,usr); 方位角[rad]
 *    posxyz sat; 衛星の位置
 *    posxyz usr; ユーザ位置
 *------------------------------------------------------------*/
double azimuth(posxyz sat,posxyz usr)
{
	posenu	enu;

	/* ENU座標に変換して方位角を求める */
	enu	=xyz_to_enu(sat,usr);
	return atan2(enu.e,enu.n);
}

/*------------------------------------------------------------
 * main() - メイン
 *------------------------------------------------------------*/
void main(int argc,char **argv)
{
	int		i;
	long	n=0L;
	double	d,data[6],data0[6],stat1[6],stat2[6],max[6],min[6];

	/* コマンドライン引数をチェック */
	if (argc!=1) {
		printf("stat - position statistics.\n");
		printf("\n");
		printf("usage: stat <in >out\n");
		printf("\n");
		exit(0);
	}

	/* 初期化 */
	for(i=0;i<6;i++) {
		stat1[i]=0.0;
		stat2[i]=0.0;
	}

	/* データファイルを読み込む */
	while(read_line(stdin)) {
		/* コメント行は読み飛ばす */
		if (linebuf[0]=='#' || linebuf[0]=='%') continue;

		d		=atof(get_field(0));			/* 時刻 */
		data[0]	=atof(get_field(0));			/* X座標[m] */
		data[1]	=atof(get_field(0));			/* y座標[m] */
		data[2]	=atof(get_field(0));			/* Z座標[m] */
		/* 水平距離 */
		data[3]	=sqrt(SQ(data[0])+SQ(data[1]));
		/* 垂直距離 */
		data[4]	=fabs(data[2]);
		/* 幾何距離 */
		data[5]	=sqrt(SQ(data[3])+SQ(data[2]));
		/* 異常値ならスキップ */
		if (d==0.0 && data[5]==0.0) continue;

		/* 和の計算 */
		for(i=0;i<6;i++) {
			if (n<1) {
				data0[i]=data[i];
				data[i]	=0.0;
			} else {
				data[i]	-=data0[i];
			}
			stat1[i]+=data[i];
			stat2[i]+=SQ(data[i]);
			if (n<1 || data[i]>max[i]) max[i]=data[i];
			if (n<1 || data[i]<min[i]) min[i]=data[i];
		}
		n++;
	}

	/* データ数 */
	printf("NUM,%d\n",n);
	if (n<1) exit(0);

	/* データ数で割る */
	for(i=0;i<6;i++) {
		stat1[i]/=n;
		stat2[i]/=n;
	}

	/* 平均値 */
	printf("AVG");
	for(i=0;i<6;i++) {
		printf(",%.4f",stat1[i]+data0[i]);
	}
	printf("\n");

	/* 標準偏差 */
	printf("STD");
	for(i=0;i<6;i++) {
		printf(",%.4f",sqrt(stat2[i]-SQ(stat1[i])));
	}
	printf("\n");

	/* RMS */
	printf("RMS");
	for(i=0;i<6;i++) {
		printf(",%.4f",
			sqrt(SQ(data0[i])+2.0*data0[i]*stat1[i]+stat2[i]));
	}
	printf("\n");

	/* 最大値 */
	printf("MAX");
	for(i=0;i<6;i++) {
		printf(",%.4f",max[i]+data0[i]);
	}
	printf("\n");

	/* 最小値 */
	printf("MIN");
	for(i=0;i<6;i++) {
		printf(",%.4f",min[i]+data0[i]);
	}
	printf("\n");

	/* 平均位置 */
	printf("OFF");
	for(i=0;i<3;i++) {
		printf(",%.4f",stat1[i]+data0[i]);
	}
	d=SQ(stat1[0]+data0[0])+SQ(stat1[1]+data0[1]);
	printf(",%.4f",sqrt(d));
	printf(",%.4f",fabs(stat1[2]+data0[2]));
	printf(",%.4f",sqrt(d+SQ(stat1[2]+data0[2])));
	printf("\n");

	exit(0);
}

