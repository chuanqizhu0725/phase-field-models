
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <typeinfo>

using namespace std;

#define DRND(x) ((double)(x) / RAND_MAX * rand()) //乱数の設定

#define ND 200 //差分計算における計算領域一辺の分割数
#define N 3    //考慮する結晶方位の数＋１(MPF0.cppと比較して、この値を大きくしている)

int nd = ND,
    ndm = ND - 1;            //計算領域の一辺の差分分割数(差分ブロック数), ND-1を定義
int nm = N - 1, nmm = N - 2; //考慮する結晶方位の数、N-2（考慮する結晶方位の数－１）を定義
double PI = 3.141592;        //π、計算カウント数
double RR = 8.3145;          //ガス定数

double theta, theta0;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidxx, phidyy, phidxy;
double ep, ep1p, ep2p;

double astre;

double phi[N][ND][ND], phi2[N][ND][ND]; //フェーズフィールド、フェーズフィールド補助配列
double aij[N][N];                       //勾配エネルギー係数
double anij[N][N];
double thij[N][N];
double wij[N][N];      //ペナルティー項の係数
double mij[N][N];      //粒界の易動度
double fij[N][N];      //粒界移動の駆動力
int phiIdx[N][ND][ND]; //位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の番号
int phiNum[ND][ND];
int phinum;

int i, j, k, l, ii, jj, kk, ll, it;     //整数
int ip, ipp, im, imm, jp, jpp, jm, jmm; //整数
int n1, n2, n3;                         //整数

int istep = 0;
// int n000;		//位置(i,j)において、pが０ではない方位の個数（n00>=n000）
int nstep; //計算カウント数の最大値（計算終了カウント）
int pstep;
double dtime, L, dx;     // L計算領域の一辺の長さ(nm), 差分プロック１辺の長さ(m)
double M0;               //粒界の易動度
double W0;               //ペナルティー項の係数
double A0;               //勾配エネルギー係数
double F0;               //粒界移動の駆動力
double temp;             //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;            //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;    //モル体積

double A_alph, A_beta, c_alph0, c_beta0; //化学的自由エネルギー内のパラメータ
double c, sa, sb, Da, Db, D;             //濃度場、Sα、Sβ、拡散係数（Dα、Dβ）、重み付拡散係数D
double gii, gjj, cii, cjj;               //化学的自由エネルギー、局所平衡濃度
double ptl;
double dF;
double dgdc;                           //化学的自由エネルギーの濃度微分
double cddtt;                          //濃度の変化率
double dc0;                            //濃度場補正用の作業変数
double dev1_a, dev1_b, dev2_a, dev2_b; //拡散方程式内における微分計算の際の作業変数

double c0;                       //平均組成
double ch[ND][ND], ch2[ND][ND];  //濃度場、濃度場の補助配列
double cah[ND][ND], cbh[ND][ND]; //局所平衡（平行接線則）の濃度場
double sah[ND][ND], sbh[ND][ND]; //α相とβ相の存在確率（SαとSβ）

void initialize();
void datasave(int step);
void datain();
void log();
double calcTheta(double dy, double dx);

//************ 初期場(フェーズフィールド)の設定サブルーチン *************
void initialize()
{
    int i, j, k, l, it;           //整数
    int ii, jj, kk;               //整数
    int ip, im, jp, jm;           //整数
    int x1, y1, x1h[10], y1h[10]; //初期核の座標
    double sum1, t, r0, r;
    srand(3.0); // 乱数初期化
                // srand(time(NULL)); // 乱数初期化

    x1h[1] = 0.25 * nd;
    y1h[1] = 0.25 * nd; //初期核１の座標設定
    // x1h[2] = 0.75 * nd;
    // y1h[2] = 0.75 * nd; //初期核２の座標設定
    // x1h[1] = 0.25 * nd;
    // y1h[1] = 0.25 * nd; //初期核３の座標設定

    //*** 式(4.36) - 式(4.39)の配列（K,W,M,E）の設定 **************************
    for (i = 1; i <= nm; i++)
    {
        for (j = 1; j <= nm; j++)
        {
            wij[i][j] = W0;
            aij[i][j] = A0;
            mij[i][j] = M0;
            fij[i][j] = 0.0;
            anij[i][j] = false;
            thij[i][j] = 0.0;
            if ((i == nm) || (j == nm))
            {
                fij[i][j] = F0;
                anij[i][j] = true;
                thij[i][j] = PI / 2.0 * DRND(1);
            }
            if (i > j)
            {
                fij[i][j] = -fij[i][j];
                thij[i][j] = thij[j][i];
            }
            if (i == j)
            {
                wij[i][j] = 0.0;
                aij[i][j] = 0.0;
                mij[i][j] = 0.0;
                fij[i][j] = 0.0;
                thij[i][j] = 0.0;
                anij[i][j] = false;
            }
        }
    }

    //*** 初期場の設定 *****************************************
    for (k = 1; k <= nm; k++)
    {
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                phi[k][i][j] = 0.0;
                phi2[k][i][j] = 0.0; //フェーズフィールド、および補助配列の初期化
            }
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            for (ii = 1; ii <= nm - 1; ii++)
            {
                phi[ii][i][j] = 0.0;
            }
            phi[nm][i][j] = 1.0; // nm番目のフェーズフィールドを１に初期化
        }
    }

    r0 = 15.0;
    for (ii = 1; ii <= nm - 1; ii++)
    {
        x1 = x1h[ii];
        y1 = y1h[ii];
        // x1 = nd * DRND(1);
        // y1 = nd * DRND(1); //初期核の位置
        for (i = 0; i <= ndm; i++)
        {
            for (j = 0; j <= ndm; j++)
            {
                r = sqrt((double(i - x1)) * (double(i - x1)) + (double(j - y1)) * (double(j - y1)));
                if (r <= r0)
                {
                    phi[ii][i][j] = 1.0;
                    phi[nm][i][j] = 0.0;
                } //初期核位置のフェーズフィールドを設定
            }
        }
    }

    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            sum1 = 0.0;
            for (k = 1; k <= nm; k++)
            {
                sum1 += phi[k][i][j];
            }
            for (k = 1; k <= nm; k++)
            {
                phi[k][i][j] = phi[k][i][j] / sum1;
            } //フェーズフィールドの規格化補正
        }
    }

    double count_l = 0.0;
    double count_s = 0.0;
    for (i = 0; i <= ndm; i++)
    {
        for (j = 0; j <= ndm; j++)
        {
            if (phi[nm][i][j] == 1.0)
            {
                count_l += 1.0;
                ch[i][j] = 0.2;
            }
            else
            {
                count_s += 1.0;
                ch[i][j] = 0.1;
            }
        }
    }

    c0 = count_l / (count_s + count_l) * 0.2 + count_s / (count_s + count_l) * 0.1;
    cout << "count_l: " << count_l << endl;
    cout << "count_s: " << count_s << endl;
    cout << "nominal concentration is: " << c0 << endl;
}

void datasave(int step)
{
    cout << "hello, I started one step!" << step << "\n";

    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        for (int j = 0; j <= ndm; j++)
        {
            double sumPic = 0.0;
            for (int k = 0; k <= nm; k++)
            {
                sumPic += phi[k][i][j] * phi[k][i][j];
            }
            fprintf(stream, "%e   ", sumPic); //フェーズフィールドの保存
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ

    FILE *stream1; //ストリームのポインタ設定
    char buffer1[30];
    sprintf(buffer1, "data/con/2d%d.csv", step);
    stream1 = fopen(buffer1, "a"); //書き込む先のファイルを追記方式でオープン

    for (int i = 0; i <= ndm; i++)
    {
        for (int j = 0; j <= ndm; j++)
        {
            fprintf(stream1, "%e   ", ch[i][j]); //フェーズフィールドの保存
            fprintf(stream1, "\n");
        }
    }
    fclose(stream1); //ファイルをクローズ
}

void datain()
{
    // Create a text string, which is used to output the text file
    string lineText;
    string paraText;
    string dataText;
    string delimiter = "#";

    // Read from the text file
    ifstream inputfile("input.txt");

    // Use a while loop together with the getline() function to read the file line by line
    while (getline(inputfile, lineText))
    {
        // Output the text from the file
        paraText = lineText.substr(0, lineText.find(delimiter));
        dataText = lineText.substr(lineText.find(delimiter) + 1);
        if (paraText == "nstep")
        {
            nstep = stoi(dataText);
        }
        else if (paraText == "pstep")
        {
            pstep = stoi(dataText);
        }
        else if (paraText == "dtime")
        {
            dtime = stod(dataText);
        }
        else if (paraText == "temp")
        {
            temp = stod(dataText);
        }
        else if (paraText == "L")
        {
            L = stod(dataText);
        }
        else if (paraText == "vm0")
        {
            vm0 = stod(dataText);
        }
        else if (paraText == "delta")
        {
            delta = stod(dataText);
        }
        else if (paraText == "mobi")
        {
            mobi = stod(dataText);
        }
        else if (paraText == "astre")
        {
            astre = stod(dataText);
        }
    }
    // Close the file
    inputfile.close();
}

void log()
{
    cout << "-----------------------\n"
         << "grid length: " << dx << " m\n"
         << "interface energy: " << gamma0 << " non-dimen unit\n"
         << "anisotropy strength: " << astre << "\n";
}

double calcTheta(double dy, double dx)
{
    if (dx != 0.0)
    {
        return atan(dy / dx);
    }
    else if (dx == 0.0 && dy > 0)
    {
        return PI / 2.0;
    }
    else if (dx == 0.0 && dy < 0)
    {
        return PI / (-2.0);
    }
    else
    {
        return 0;
    }
}