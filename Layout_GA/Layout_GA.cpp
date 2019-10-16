#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstring>
#include <algorithm>

#include <opencv2/opencv.hpp>

#include "Layout_GA.h"
 
using namespace std;
   
#define POPSIZE 50//种群内个体数量
#define MAXGENS 999//最大的迭代次数

//交叉概率和突变概率有自适应改变方法, 即后面的dCrossP 和 dMutateP
//#define PXOVER 0.8//交叉概率
//#define PMUTATION 0.15//突变概率
#define PXOVER 0.8//交叉概率
#define PMUTATION 0.2//突变概率

#define BIG_BOARD_WIDTH 2000//大板宽
#define BIG_BOARD_HEIGHT 1000//大板高
#define MARGIN 6 //在实际生产中很有意义，比如切割玻璃或皮革时要把边角料切掉
#define GAP 3  //切割的毛坯之间的间隔

INDIVIDUAL tPopulation[POPSIZE + 1];//最后一个用于保存适应度最大的个体
INDIVIDUAL tNewPopulation[POPSIZE + 1];

MYRECT tSmallRects[111];//存储要切割的矩形毛坯信息
double dSmallRects[112][2];//存储要切割矩形排序之前的信息，从下标1开始保存，第一列为矩形宽

int smallRectsNum;//要切割的个数，即读入txt的数据行数
 
int main();
void crossover();//交叉操作
void elitist();//精英策略，即保留适应度最高的
void evaluate();//适应度函数，在这里是排样最大高度的倒数
int Int_uniform_ab(int a, int b);//[a,b]之间生成随机整数
void initialize(string filename);//初始化种群
void keep_the_best();//保存最优，即保存适应度最高
void mutate();//变异
double Dou_uniform_ab(double a, double b);//[a,b]之间生成随机double型数据
void report(int generation);//输出算法进度
void selector();//选代
void timestamp();//在控制台输出中打印时间戳，于遗传算法无实际意义
void Xover(int one, int two);//对选择的两个个体进行单点交叉
void Xover2(int one, int two);//对选择的两个个体进行双点交叉
double get_variation_coefficient();//该系数可以表示种群数据集中程度，用于自适应计算交叉概率Pc和变异概率Pm
 
/*********************************************************************************************************
流程：在每次迭代中，选择适应度高的进行遗传，这个过程是通过轮盘赌算法来实现，
即适应度高的占的比例大，选中的概率相应也大，然后进行交叉和变异操作，再根据
适应度函数进行评价，直到迭代完毕
*********************************************************************************************************/

int main()
{
	srand((unsigned)time(NULL));
 
	timestamp();

	cout << "\nLAYOYT_GA: \n" << "  C++ version \n";
	cout << "  A layout example of GA.\n";
 
	string filename = "./data/layout_ga_input.txt";//要切割的矩形毛坯的宽和高
	initialize(filename);//在txt中读出要切割的尺寸信息，适应值全部初始为0
 
	evaluate();//适应度函数
 
	keep_the_best();//找出适应度值最高的个体，存放在种群数组末尾
 
	for (int generation = 0; generation < MAXGENS; ++generation )
	{
		selector();//轮盘赌算法
		crossover();//单点交叉，即把交叉点之前换掉
		mutate();//
		report(generation);
		evaluate();
		elitist();
	}
 
	cout << "\n  Best member after " << MAXGENS << " generations:\n\n";
 
	for (int i = 0; i < smallRectsNum; ++i )
	{
		if(0 == (i + 1) % 4)
		{
			cout << setw(11) << "var("<< i <<"): "<< tPopulation[POPSIZE].nGene[i] <<"" << endl;
		}
		else
		{
			cout << setw(11) << "var("<< i <<"): "<< tPopulation[POPSIZE].nGene[i] <<"";
		}
	}
 
	cout << "\n\n  Best fitness = " << tPopulation[POPSIZE].fitness;
	cout << "\n\nLAYOUT_GA: Normal end of execution.\n\n";

	timestamp();

	//cout << endl;
	//cout << "Best fitness --> Coordinates: " << endl;
	//for(int i = 0; i < smallRectsNum; ++i)
	//{
	//	cout << tPopulation[POPSIZE].tPosition[i].dX;
	//	cout << "  " << tPopulation[POPSIZE].tPosition[i].dY << endl;
	//}

	//绘图保存
    cv::Mat img(BIG_BOARD_HEIGHT + 222, BIG_BOARD_WIDTH + 222, CV_8UC3, cv::Scalar(0));
 
    cv::Point root_points[1][4];
	//画大板
	//Point(x, y)   x-->img.col
	root_points[0][0] = cv::Point(111, 111);
    root_points[0][1] = cv::Point(111 + BIG_BOARD_WIDTH, 111);
    root_points[0][2] = cv::Point(111 + BIG_BOARD_WIDTH, 111 + BIG_BOARD_HEIGHT);
    root_points[0][3] = cv::Point(111, 111 + BIG_BOARD_HEIGHT);

	const cv::Point* ppt[1] = {root_points[0]};
    int npt[] = {4};
    polylines(img, ppt, npt, 1, 1, cv::Scalar(255, 255, 255), 3, 8, 0);

#if 1
	//画小矩形并填充
	for(int i = 0; i < smallRectsNum; ++i)
	{
		int n = tPopulation[POPSIZE].nGene[i];
		if(n > 0)
		{
			root_points[0][0] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + 0.5) + 111);
			root_points[0][1] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + dSmallRects[n][0] + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + 0.5) + 111);
			root_points[0][2] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + dSmallRects[n][0] + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + dSmallRects[n][1] + 0.5) + 111);
			root_points[0][3] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + dSmallRects[n][1] + 0.5) + 111);
		}
		else
		{
			root_points[0][0] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + 0.5) + 111);
			root_points[0][1] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + dSmallRects[-n][1] + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + 0.5) + 111);
			root_points[0][2] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + dSmallRects[-n][1] + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + dSmallRects[-n][0] + 0.5) + 111);
			root_points[0][3] = cv::Point((int)(tPopulation[POPSIZE].tPosition[i].dY + 0.5) + 111,
				(int)(tPopulation[POPSIZE].tPosition[i].dX + dSmallRects[-n][0] + 0.5) + 111);
		}

		const cv::Point* ppt[1] = {root_points[0]};
		int npt[] = {4};
		polylines(img, ppt, npt, 1, 1, cv::Scalar(255, 255, 255), 1, 8, 0);
		fillPoly(img, ppt, npt, 1, cv::Scalar(Int_uniform_ab(0, 255),
			Int_uniform_ab(0, 255), Int_uniform_ab(0, 255)));
	}

#endif
    imshow("LAYOUT_GA", img);
    cv::waitKey();
	cv::imwrite("GA_RES.jpg", img);

	system("pause");
 
	return 0;
}


double get_variation_coefficient()
{
#if 0
	double dFitnessSum = 0.0;
	double dFitnessQuadraticSum = 0.0;
	for(int i = 0; i < POPSIZE; ++i)
	{
		dFitnessSum += tPopulation[i].fitness;
		dFitnessQuadraticSum += tPopulation[i].fitness * tPopulation[i].fitness;
	}
	double EX = dFitnessSum / (double)POPSIZE;
	double DX = dFitnessQuadraticSum / (double)POPSIZE - EX * EX;
	
	//return (EX + 1) / sqrt(DX);
	return EX / sqrt(DX);
#endif

	double dFitnessSum = 0.0;
	for (int i = 0; i < POPSIZE; ++i)
	{
		dFitnessSum += tPopulation[i].fitness;
	}
 
	double avg = dFitnessSum / (double)POPSIZE;

	double dSquareSum = 0.0;
	for(int i = 0; i < POPSIZE; ++i)
	{
		dSquareSum += (tPopulation[i].fitness - avg) * (tPopulation[i].fitness - avg);
	}
	double stddev = sqrt(dSquareSum / (POPSIZE - 1));

	return stddev / avg;
}

void crossover()//单点交叉结果貌似好一点，最好适应度0.000682594
{
	const double a = 0.0;
	const double b = 1.0;
	int one;
	int first = 0;

	//double w = get_variation_coefficient();
	//const double k = 9.99;// 0 <= k <= +无穷
	//double dCrossP = 1.0 / (1 + (exp(-w * k))) - 0.1;//逻辑函数

	bool cross_flag[POPSIZE];//记录交叉位置，交叉过的就跳过
	memset(cross_flag, false, sizeof(cross_flag));

	double f_max = 0.0;//记录种群最大适应值
	double f_sum = 0.0;//记录种群适应度总和，用以计算种群适应度平均值
	for(int i = 0; i < POPSIZE; ++i)
	{
		f_max = max(f_max, tPopulation[i].fitness);
		f_sum += tPopulation[i].fitness;
	}
	double f_avg = f_sum / (double)(POPSIZE);
		
	double p_c;//交叉概率
	//k1 < k2
	const double k1 = 0.9;
	const double k2 = 1.0;

 
	for( int mem = 0; mem < POPSIZE; ++mem )
	{
#if 1
		double x = Dou_uniform_ab(a, b);
		if(x < PXOVER)//0.8的概率交叉
		//if(x < dCrossP)
		{
			++first;
 
			if(0 == first % 2)
			{
				//Xover2(one, mem);
				Xover(one, mem);
			}
			else
			{
				one = mem;
			}
		}
#endif

#if 0
		//随机生成两个要交叉的个体下标
		int num1 = Int_uniform_ab(0, POPSIZE - 1);
		int num2 = Int_uniform_ab(0, POPSIZE - 1);
		double f = max(tPopulation[num1].fitness, tPopulation[num2].fitness);//要交叉的两个个体中较大的适应度值
		if(f >= f_avg)
		{
			p_c = k1 * (f_max - f) / (f_max - f_avg);
		}
		else
		{
			p_c = k2;
		}
		double x = Dou_uniform_ab(a, b);
		if(x < p_c && false == cross_flag[num1] && false == cross_flag[num2])
		{
			Xover(num1, num2);
			cross_flag[num1] = cross_flag[num2] = true;
		}
		#endif

	}
	return;
}
 
void elitist()//保存前代最好的，即适应度值最大的
{
	double best;
	int best_mem;
	double worst;
	int worst_mem;
 
	best = worst = tPopulation[0].fitness;
 
	for (int i = 0; i < POPSIZE - 1; ++i)
	{
		if (tPopulation[i + 1].fitness < tPopulation[i].fitness)
		{
			if (best <= tPopulation[i].fitness)
			{
				best = tPopulation[i].fitness;
				best_mem = i;
			}
 
			if (tPopulation[i + 1].fitness <= worst)
			{
				worst = tPopulation[i + 1].fitness;
				worst_mem = i + 1;
			}
		}
		else
		{
			if (tPopulation[i].fitness <= worst)
			{
				worst = tPopulation[i].fitness;
				worst_mem = i;
			}
 
			if (best <= tPopulation[i+1].fitness)
			{
				best = tPopulation[i+1].fitness;
				best_mem = i + 1;
			}
		}
	}

	if (tPopulation[POPSIZE].fitness <= best)
	{
		tPopulation[POPSIZE] = tPopulation[best_mem];
	}
	else
	{
		tPopulation[worst_mem] = tPopulation[POPSIZE];
	} 
 
	return;
}
 
void evaluate()//此处适应度值即最高水平线的大小的倒数，默认大板充足，即没有排不下的情况
{
	//dYTmp每一层更新一次：更新规则即上一层最高值 + gap
	for(int member = 0; member < POPSIZE; ++member)
	{
		double dXTmp = MARGIN;
		double dYTmp = MARGIN;
		double dHighestTmp = MARGIN;
		for(int i = 0; i < smallRectsNum; ++i)
		{
			int n = tPopulation[member].nGene[i];
			if(n > 0 && dXTmp + dSmallRects[n][1] <= BIG_BOARD_HEIGHT - MARGIN)
			{
				tPopulation[member].tPosition[i].dX = dXTmp;
				tPopulation[member].tPosition[i].dY = dYTmp;
				dXTmp += dSmallRects[n][1];
				dHighestTmp = max(dHighestTmp, dYTmp + dSmallRects[n][0]);
				dXTmp += GAP;
			}
			else if(n < 0 && dXTmp + dSmallRects[-n][0] <= BIG_BOARD_HEIGHT - MARGIN)
			{
				tPopulation[member].tPosition[i].dX = dXTmp;
				tPopulation[member].tPosition[i].dY = dYTmp;
				dXTmp += dSmallRects[-n][0];
				dHighestTmp = max(dHighestTmp, dYTmp + dSmallRects[-n][1]);
				dXTmp += GAP;
			}
			else
			{
				dXTmp = MARGIN;
				dYTmp = dHighestTmp + GAP;
				--i;
			}
		}
		tPopulation[member].fitness = 1.0 / dHighestTmp;
	}
	return;
}

int Int_uniform_ab(int a, int b)//[a, b]
{	
	return (rand() % (b - a + 1)) + a;
}

 //从文档中读取切割信息，生成初始种群
void initialize (string filename)
{
	ifstream input;
 
	input.open (filename.c_str());
 
	if(!input)
	{
		cerr << "\n";
		cerr << "INITIALIZE - Fatal error!\n";
		cerr << "  Cannot open the input file!\n";
		exit(1);
	}

	int cnt = 0;
	memset(dSmallRects, 0, sizeof(dSmallRects));
	while(!input.eof())
	{
		double dW, dH;
		input >> dW >> dH;
		tSmallRects[cnt] = MYRECT(dW, dH, cnt + 1);
		dSmallRects[cnt + 1][0] = tSmallRects[cnt].dWidth;
		dSmallRects[cnt + 1][1] = tSmallRects[cnt].dHeight;
		cnt++;
	}
	input.close();
	smallRectsNum = cnt;
	//本来sort函数默认升序，但在结构体中已重载小于号，所以现在是降序，即面积大的优先排列
	sort(tSmallRects, tSmallRects + smallRectsNum);

	memset(tPopulation, 0, sizeof(INDIVIDUAL) * (POPSIZE+ 1));
	memset(tNewPopulation, 0, sizeof(INDIVIDUAL) * (POPSIZE+ 1));

	for(int i = 0; i < POPSIZE; ++i)
	{
		for(int j = 0; j < smallRectsNum; ++j)
		{
			tPopulation[i].nGene[j] = (rand() % 2) ? (tSmallRects[j].nIndex) : (-tSmallRects[j].nIndex);
		}
	}
 
	return;
}
 
void keep_the_best()
{ 
	for ( int mem = 0; mem < POPSIZE; ++mem )
	{
		if ( tPopulation[POPSIZE].fitness < tPopulation[mem].fitness )
		{
			tPopulation[POPSIZE] = tPopulation[mem];
		}
	}
 
	return;
}

//变异策略：生成n1  n2  n3三个随机数，其中n1和n2位置互换，n3取反
void mutate()
{
	const double a = 0.0;
	const double b = 1.0;

	//double w = get_variation_coefficient();
	//cout << "w: " << w << endl;
	//const double k = 0.8;//0 <= k <= 1
	//double dMutateP = k / (5 * (1 + (exp(w))));

	double f_max = 0.0;//记录种群最大适应值
	double f_sum = 0.0;//记录种群适应度总和，用以计算种群适应度平均值
	for(int i = 0; i < POPSIZE; ++i)
	{
		f_max = max(f_max, tPopulation[i].fitness);
		f_sum += tPopulation[i].fitness;
	}
	double f_avg = f_sum / (double)(POPSIZE);
		
	double p_m;//变异概率
	//k1 < k2
	const double k1 = 0.4;
	const double k2 = 0.5;

	for(int i = 0; i < POPSIZE; ++i)
	{
		if(tPopulation[i].fitness >= f_avg)
		{
			p_m = k1 * (f_max - tPopulation[i].fitness) / (f_max - f_avg);
		}
		else
		{
			p_m = k2;
		}
		double x = Dou_uniform_ab(a, b);	 
		if(x < PMUTATION)
		//if(x < dMutateP)
		//if(x < p_m)
		{
			int n1 = Int_uniform_ab(0, smallRectsNum - 1);
			int n2 = Int_uniform_ab(0, smallRectsNum - 1);
			int t = tPopulation[i].nGene[n1];
			tPopulation[i].nGene[n1] = tPopulation[i].nGene[n2];
			tPopulation[i].nGene[n2] = t;
			int n3 = Int_uniform_ab(0, smallRectsNum - 1);
			tPopulation[i].nGene[n3] = -tPopulation[i].nGene[n3];
		}
	}
 
	return;
}

double Dou_uniform_ab(double a, double b)
{
	//rand() / double(RAND_MAX)可以生成0~1之间的浮点数
	return a + static_cast<double>(rand()) / RAND_MAX * (b - a);
}

void report(int generation)
{
	if (0 == generation)
	{
		cout << "\n";
		cout << "  Generation       Best            Average       Standard \n";
		cout << "  number           value           fitness       deviation \n";
		cout << "\n";
	}

 	double dFitnessSum = 0.0;
	for (int i = 0; i < POPSIZE; ++i )
	{
		dFitnessSum += tPopulation[i].fitness;
	}
 
	double avg = dFitnessSum / (double)POPSIZE;

	double dSquareSum = 0.0;
	for(int i = 0; i < POPSIZE; ++i)
	{
		dSquareSum += (tPopulation[i].fitness - avg) * (tPopulation[i].fitness - avg);
	}
	double stddev = sqrt(dSquareSum / (POPSIZE - 1));
	double best_val = tPopulation[POPSIZE].fitness;
 
	cout << "  " << setw(8) << generation 
		<< "  " << setw(14) << best_val 
		<< "  " << setw(14) << avg 
		<< "  " << setw(14) << stddev << "\n";
 
	return;
}
 
void selector()
{
	const double a = 0.0;
	const double b = 1.0;
	double sum = 0.0;

	for ( int mem = 0; mem < POPSIZE; ++mem )
	{
		sum = sum + tPopulation[mem].fitness;
	}

	for ( int mem = 0; mem < POPSIZE; ++mem )
	{
		tPopulation[mem].rfitness = tPopulation[mem].fitness / sum;
	}
	
	tPopulation[0].cfitness = tPopulation[0].rfitness;
	for ( int mem = 1; mem < POPSIZE; ++mem )
	{
		tPopulation[mem].cfitness = tPopulation[mem - 1].cfitness +       
			tPopulation[mem].rfitness;
	}

	for ( int i = 0; i < POPSIZE; ++i )
	{ 
		double p = Dou_uniform_ab(a, b);//产生（0,1）区间一个double型数据
		if ( p < tPopulation[0].cfitness )
		{
			tNewPopulation[i] = tPopulation[0];      
		}
		else
		{
			for ( int j = 0; j < POPSIZE; ++j )
			{ 
				if ( tPopulation[j].cfitness <= p && p < tPopulation[j + 1].cfitness )
				{
					tNewPopulation[i] = tPopulation[j + 1];
				}
			}
		}
	}
	
	for ( int i = 0; i < POPSIZE; ++i )
	{
		tPopulation[i] = tNewPopulation[i]; 
	}

	return;     
}

void timestamp()
{
# define TIME_SIZE 40
 
	static char time_buffer[TIME_SIZE];
	const struct tm *tm;
	size_t len;
	time_t now;
 
	now = time (NULL);
	tm = localtime (&now);
	//将时间格式化
	len = strftime (time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
 
	cout << time_buffer << "\n";
 
	return;
# undef TIME_SIZE
}
 
void Xover (int one, int two)
{
	int nPoint = Int_uniform_ab(0, smallRectsNum - 1);

	//交叉后的基因是不合法的
	int * pInt2 = new int [nPoint];//交叉操作后，存放基因2重复的序列
	int * pInt1 = new int [nPoint];//交叉操作后，存放基因1重复的序列

	//memset(pInt1, 0, nPoint * sizeof(pInt1));
	//memset(pInt2, 0, nPoint * sizeof(pInt2));

	int j;
	int cnt = 0;
	for(int i = 0; i < nPoint; ++i)
	{
		for(j = 0; j < nPoint; ++j)
		{
			if(abs(tPopulation[one].nGene[i]) == abs(tPopulation[two].nGene[j]))
			{
				break;
			}
		}
		if(nPoint == j)
		{
			pInt2[cnt++] = tPopulation[one].nGene[i];
		}
	}

	cnt = 0;
	for(int i = 0; i < nPoint; ++i)
	{
		for(j = 0; j < nPoint; ++j)
		{
			if(abs(tPopulation[one].nGene[j]) == abs(tPopulation[two].nGene[i]))
			{
				break;
			}
		}
		if(nPoint == j)
		{
			pInt1[cnt++] = tPopulation[two].nGene[i];
		}
	}

	//把基因2中重复序列的所在位置记录下来，依次赋值给基因1重复序列处
	int * pIntRep2 = new int [cnt];//存放基因2重复的序列，按在基因2中出现顺序
	int * pIntRep1 = new int [cnt];//存放基因1重复的序列，按在基因1中出现顺序

	int num = 0;

	for(int i = 0; i < smallRectsNum; ++i)
	{
		for(j = 0; j < cnt; ++j)
		{
			if(abs(pInt2[j]) == abs(tPopulation[two].nGene[i]))
			{
				pIntRep2[num++] = tPopulation[two].nGene[i];
			}
		}
		if(cnt == num)//
		{
			break;
		}
	}

	num = 0;
	for(int i = 0; i < smallRectsNum; ++i)
	{
		for(j = 0; j < cnt; ++j)
		{
			if(abs(pInt1[j]) == abs(tPopulation[one].nGene[i]))
			{
				pIntRep1[num++] = tPopulation[one].nGene[i];
			}
		}
		if(cnt == num)
		{
			break;
		}
	}

		//修改不合法基因1
	for(int i = 0; i < cnt; ++i)
	{
		for(j = nPoint; j < smallRectsNum; ++j)
		{
			if(tPopulation[one].nGene[j] == pIntRep1[i])
			{
				tPopulation[one].nGene[j] = pIntRep2[i];
				break;
			}
		}
	}

	//修改不合法基因2
	for(int i = 0; i < cnt; ++i)
	{
		for(j = nPoint; j < smallRectsNum; ++j)
		{
			if(tPopulation[two].nGene[j] == pIntRep2[i])
			{
				tPopulation[two].nGene[j] = pIntRep1[i];
				break;
			}
		}
	}

	for ( int i = 0; i < nPoint; ++i )
	{
		int t = tPopulation[one].nGene[i];
		tPopulation[one].nGene[i] = tPopulation[two].nGene[i];
		tPopulation[two].nGene[i] = t;
	}

	if(pInt2 != NULL)
	{
		delete [] pInt2;
	}
	if(pInt1 != NULL)
	{
		delete [] pInt1;
	}
	if(pIntRep2 != NULL)
	{
		delete [] pIntRep2;
	}
	if(pIntRep1 != NULL)
	{
		delete [] pIntRep1;
	}

	return;
}

void Xover2 (int one, int two)//两点交叉
{
	int nPoint1 = Int_uniform_ab(0, smallRectsNum - 1);
	int nPoint2 = Int_uniform_ab(0, smallRectsNum - 1);

	if(nPoint1 > nPoint2)
	{
		int t = nPoint1;
		nPoint1 = nPoint2;
		nPoint2 = t;
	}

	//交叉后的基因是不合法的
	int * pInt2 = new int [nPoint2 - nPoint1];//交叉操作后，存放基因2重复的序列
	int * pInt1 = new int [nPoint2 - nPoint1];//交叉操作后，存放基因1重复的序列

	int j;
	int cnt = 0;
	for(int i = nPoint1; i < nPoint2; ++i)
	{
		for(j = nPoint1; j < nPoint2; ++j)
		{
			if(abs(tPopulation[one].nGene[i]) == abs(tPopulation[two].nGene[j]))
			{
				break;
			}
		}
		if(nPoint2 == j)
		{
			pInt2[cnt++] = tPopulation[one].nGene[i];
		}
	}

	cnt = 0;
	for(int i = nPoint1; i < nPoint2; ++i)
	{
		for(j = nPoint1; j < nPoint2; ++j)
		{
			if(abs(tPopulation[one].nGene[j]) == abs(tPopulation[two].nGene[i]))
			{
				break;
			}
		}
		if(nPoint2 == j)
		{
			pInt1[cnt++] = tPopulation[two].nGene[i];
		}
	}

	//把基因2中重复序列的所在位置记录下来，依次赋值给基因1重复序列处
	int * pIntRep2 = new int [cnt];//存放基因2重复的序列，按在基因2中出现顺序
	int * pIntRep1 = new int [cnt];//存放基因1重复的序列，按在基因1中出现顺序

	int num = 0;

	for(int i = 0; i < smallRectsNum; ++i)
	{
		for(j = 0; j < cnt; ++j)
		{
			if(abs(pInt2[j]) == abs(tPopulation[two].nGene[i]))
			{
				pIntRep2[num++] = tPopulation[two].nGene[i];
			}
		}
		if(cnt == num)//
		{
			break;
		}
	}

	num = 0;
	for(int i = 0; i < smallRectsNum; ++i)
	{
		for(j = 0; j < cnt; ++j)
		{
			if(abs(pInt1[j]) == abs(tPopulation[one].nGene[i]))
			{
				pIntRep1[num++] = tPopulation[one].nGene[i];
			}
		}
		if(cnt == num)
		{
			break;
		}
	}

		//修改不合法基因1
	for(int i = 0; i < cnt; ++i)
	{
		for(j = 0; j < smallRectsNum; ++j)
		{
			if(tPopulation[one].nGene[j] == pIntRep1[i])
			{
				tPopulation[one].nGene[j] = pIntRep2[i];
				break;
			}
		}
	}

	//修改不合法基因2
	for(int i = 0; i < cnt; ++i)
	{
		for(j = 0; j < smallRectsNum; ++j)
		{
			if(tPopulation[two].nGene[j] == pIntRep2[i])
			{
				tPopulation[two].nGene[j] = pIntRep1[i];
				break;
			}
		}
	}

	for ( int i = nPoint1; i < nPoint2; ++i )
	{
		int t = tPopulation[one].nGene[i];
		tPopulation[one].nGene[i] = tPopulation[two].nGene[i];
		tPopulation[two].nGene[i] = t;
	}

	if(pInt2 != NULL)
	{
		delete [] pInt2;
	}
	if(pInt1 != NULL)
	{
		delete [] pInt1;
	}
	if(pIntRep2 != NULL)
	{
		delete [] pIntRep2;
	}
	if(pIntRep1 != NULL)
	{
		delete [] pIntRep1;
	}

	return;
}
