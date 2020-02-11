/*------------------------------------------------------------
 * TEST1.c - Practice for Position Computation.
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
 * 測位計算
 *------------------------------------------------------------*/

/*------------------------------------------------------------
 * inverse_matrix() - 逆行列を計算する
 *------------------------------------------------------------
 *  void inverse_matrix(a,n);
 *    double a[][]; 元の行列(逆行列で上書きされる)
 *    int m;        行列の次元(1〜MAX_M)
 *------------------------------------------------------------
 *   与えられた行列の逆行列を求める. 結果により、元の行列が
 * 上書きされる.
 *------------------------------------------------------------*/
static void inverse_matrix(double a[MAX_M][MAX_M],int m)
{
	int		i,j,k;
	double	b[MAX_M][MAX_M+MAX_M];

	/* 操作用の行列をつくる */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			b[i][j]=a[i][j];
			if (i==j) b[i][j+m]=1.0; else b[i][j+m]=0.0;
		}
	}

	/* ガウスの消去法 */
	for(i=0;i<m;i++) {
		/* 第i行をb[i][i]で正規化する */
		if (fabs(b[i][i])<=1E-10) {
			fprintf(stderr,"Cannot inverse matrix.\n");
			exit(2);
		}
		for(j=m+m-1;j>=i;j--) {
			b[i][j]/=b[i][i];
		}

		/* 他の行の第i列を消去する */
		for(k=0;k<m;k++) if (k!=i) {
			for(j=m+m-1;j>=i;j--) {
				b[k][j]-=b[k][i]*b[i][j];
			}
		}
	}

	/* 元の行列を逆行列で上書きする */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			a[i][j]=b[i][j+m];
		}
	}
}

/*------------------------------------------------------------
 * compute_solution() - 最小二乗法で方程式を解く
 *------------------------------------------------------------
 *  void compute_solution(G,dr,wgt,dx,cov,n,m);
 *    double G[][];   デザイン行列(n×m)
 *    double dr[];    方程式の右辺(n次)
 *    double wgt[];   重み係数(n次)/NULL:重みなし
 *    double dx[];    方程式の解で上書きされる(m次)
 *    double cov[][]; 共分散行列で上書きされる(m×m)
 *    int n;          方程式の数
 *    int m;          未知数の数
 *------------------------------------------------------------
 *   与えられた方程式を最小二乗法により解く. 重みが不要な場合
 * は wgt=NULL として呼び出す.
 *------------------------------------------------------------*/
void compute_solution(double G[MAX_N][MAX_M],double dr[MAX_N],
	double wgt[MAX_N],double dx[MAX_M],double cov[MAX_M][MAX_M],
	int n,int m)
{
	int		i,j,k;
	double	w,a[MAX_M][MAX_N];

	/* GtGを求める */
	for(i=0;i<m;i++) {
		for(j=0;j<m;j++) {
			cov[i][j]=0.0;
			for(k=0;k<n;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[k];
				cov[i][j]+=G[k][i]*G[k][j]*w;
			}
		}
	}

	/* 逆行列を求める（これが共分散行列Cとなる） */
	inverse_matrix(cov,m);

	/* Gtをかける */
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			a[i][j]=0.0;
			for(k=0;k<m;k++) {
				if (wgt==NULL) w=1.0; else w=wgt[j];
				a[i][j]+=cov[i][k]*G[j][k]*w;
			}
		}
	}

	/* drをかけると解になる */
	for(i=0;i<m;i++) {
		dx[i]=0.0;
		for(k=0;k<n;k++) {
			dx[i]+=a[i][k]*dr[k];
		}
	}
}

/*------------------------------------------------------------
 * main() - メイン
 *------------------------------------------------------------*/
#define	LOOP	8
#define	SATS	5
static posxyz	position[SATS]={
	{-13897607.6294,-10930188.6233,19676689.6804},	/* PRN 05 */
	{-17800899.1998,15689920.8120,11943543.3888},	/* PRN 14 */
	{-1510958.2282,26280096.7818,-3117646.1949},	/* PRN 16 */
	{-12210758.3517,20413597.0201,-11649499.5474},	/* PRN 22 */
	{-170032.6981,17261822.6784,20555984.4061},		/* PRN 25 */
};
static double	range[SATS]={
	23634878.5219,	/* PRN 05 */
	20292688.3557,	/* PRN 14 */
	24032055.0372,	/* PRN 16 */
	24383229.3740,	/* PRN 22 */
	22170992.8178,	/* PRN 25 */
};

void main(int argc,char **argv)
{
	int		i,n,loop;
	double	r,G[MAX_N][MAX_M],dr[MAX_N],dx[MAX_M];
	double	sol[MAX_M],cov[MAX_M][MAX_M];
	posxyz	satpos;

	/* 解を初期化 */
	for(i=0;i<MAX_M;i++) sol[i]=0.0;

	/* 解を求めるループ */
	for(loop=0;loop<LOOP;loop++) {
		n=SATS;
		for(i=0;i<n;i++) {
			satpos	=position[i];

			/* デザイン行列をつくる */
			r		=sqrt((satpos.x-sol[0])*(satpos.x-sol[0])
						+(satpos.y-sol[1])*(satpos.y-sol[1])
						+(satpos.z-sol[2])*(satpos.z-sol[2]));
			G[i][0]	=(sol[0]-satpos.x)/r;
			G[i][1]	=(sol[1]-satpos.y)/r;
			G[i][2]	=(sol[2]-satpos.z)/r;

			/* 擬似距離の修正量 */
			dr[i]	=range[i]-r;
		}

		/* 方程式を解く */
		compute_solution(G,dr,NULL,dx,cov,n,3);

		/* 初期値に加える */
		for(i=0;i<3;i++) {
			sol[i]+=dx[i];
		}

		/* 途中経過を出力する */
		printf("LOOP %d: X=%.3f, Y=%.3f, Z=%.3f\n",
			loop+1,sol[0],sol[1],sol[2]);
	}

	exit(0);
}

