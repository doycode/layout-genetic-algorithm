#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstring>
#include <algorithm>

#include <opencv2/opencv.hpp>

#include "Layout_GA.h"
 
using namespace std;
   
#define POPSIZE 50//��Ⱥ�ڸ�������
#define MAXGENS 999//���ĵ�������

//������ʺ�ͻ�����������Ӧ�ı䷽��, �������dCrossP �� dMutateP
//#define PXOVER 0.8//�������
//#define PMUTATION 0.15//ͻ�����
#define PXOVER 0.8//�������
#define PMUTATION 0.2//ͻ�����

#define BIG_BOARD_WIDTH 2000//����
#define BIG_BOARD_HEIGHT 1000//����
#define MARGIN 6 //��ʵ�������к������壬�����и����Ƥ��ʱҪ�ѱ߽����е�
#define GAP 3  //�и��ë��֮��ļ��

INDIVIDUAL tPopulation[POPSIZE + 1];//���һ�����ڱ�����Ӧ�����ĸ���
INDIVIDUAL tNewPopulation[POPSIZE + 1];

MYRECT tSmallRects[111];//�洢Ҫ�и�ľ���ë����Ϣ
double dSmallRects[112][2];//�洢Ҫ�и��������֮ǰ����Ϣ�����±�1��ʼ���棬��һ��Ϊ���ο�

int smallRectsNum;//Ҫ�и�ĸ�����������txt����������
 
int main();
void crossover();//�������
void elitist();//��Ӣ���ԣ���������Ӧ����ߵ�
void evaluate();//��Ӧ�Ⱥ��������������������߶ȵĵ���
int Int_uniform_ab(int a, int b);//[a,b]֮�������������
void initialize(string filename);//��ʼ����Ⱥ
void keep_the_best();//�������ţ���������Ӧ�����
void mutate();//����
double Dou_uniform_ab(double a, double b);//[a,b]֮���������double������
void report(int generation);//����㷨����
void selector();//ѡ��
void timestamp();//�ڿ���̨����д�ӡʱ��������Ŵ��㷨��ʵ������
void Xover(int one, int two);//��ѡ�������������е��㽻��
void Xover2(int one, int two);//��ѡ��������������˫�㽻��
double get_variation_coefficient();//��ϵ�����Ա�ʾ��Ⱥ���ݼ��г̶ȣ���������Ӧ���㽻�����Pc�ͱ������Pm
 
/*********************************************************************************************************
���̣���ÿ�ε����У�ѡ����Ӧ�ȸߵĽ����Ŵ������������ͨ�����̶��㷨��ʵ�֣�
����Ӧ�ȸߵ�ռ�ı�����ѡ�еĸ�����ӦҲ��Ȼ����н���ͱ���������ٸ���
��Ӧ�Ⱥ����������ۣ�ֱ���������
*********************************************************************************************************/

int main()
{
	srand((unsigned)time(NULL));
 
	timestamp();

	cout << "\nLAYOYT_GA: \n" << "  C++ version \n";
	cout << "  A layout example of GA.\n";
 
	string filename = "./data/layout_ga_input.txt";//Ҫ�и�ľ���ë���Ŀ�͸�
	initialize(filename);//��txt�ж���Ҫ�и�ĳߴ���Ϣ����Ӧֵȫ����ʼΪ0
 
	evaluate();//��Ӧ�Ⱥ���
 
	keep_the_best();//�ҳ���Ӧ��ֵ��ߵĸ��壬�������Ⱥ����ĩβ
 
	for (int generation = 0; generation < MAXGENS; ++generation )
	{
		selector();//���̶��㷨
		crossover();//���㽻�棬���ѽ����֮ǰ����
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

	//��ͼ����
    cv::Mat img(BIG_BOARD_HEIGHT + 222, BIG_BOARD_WIDTH + 222, CV_8UC3, cv::Scalar(0));
 
    cv::Point root_points[1][4];
	//�����
	//Point(x, y)   x-->img.col
	root_points[0][0] = cv::Point(111, 111);
    root_points[0][1] = cv::Point(111 + BIG_BOARD_WIDTH, 111);
    root_points[0][2] = cv::Point(111 + BIG_BOARD_WIDTH, 111 + BIG_BOARD_HEIGHT);
    root_points[0][3] = cv::Point(111, 111 + BIG_BOARD_HEIGHT);

	const cv::Point* ppt[1] = {root_points[0]};
    int npt[] = {4};
    polylines(img, ppt, npt, 1, 1, cv::Scalar(255, 255, 255), 3, 8, 0);

#if 1
	//��С���β����
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

void crossover()//���㽻����ò�ƺ�һ�㣬�����Ӧ��0.000682594
{
	const double a = 0.0;
	const double b = 1.0;
	int one;
	int first = 0;

	//double w = get_variation_coefficient();
	//const double k = 9.99;// 0 <= k <= +����
	//double dCrossP = 1.0 / (1 + (exp(-w * k))) - 0.1;//�߼�����

	bool cross_flag[POPSIZE];//��¼����λ�ã�������ľ�����
	memset(cross_flag, false, sizeof(cross_flag));

	double f_max = 0.0;//��¼��Ⱥ�����Ӧֵ
	double f_sum = 0.0;//��¼��Ⱥ��Ӧ���ܺͣ����Լ�����Ⱥ��Ӧ��ƽ��ֵ
	for(int i = 0; i < POPSIZE; ++i)
	{
		f_max = max(f_max, tPopulation[i].fitness);
		f_sum += tPopulation[i].fitness;
	}
	double f_avg = f_sum / (double)(POPSIZE);
		
	double p_c;//�������
	//k1 < k2
	const double k1 = 0.9;
	const double k2 = 1.0;

 
	for( int mem = 0; mem < POPSIZE; ++mem )
	{
#if 1
		double x = Dou_uniform_ab(a, b);
		if(x < PXOVER)//0.8�ĸ��ʽ���
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
		//�����������Ҫ����ĸ����±�
		int num1 = Int_uniform_ab(0, POPSIZE - 1);
		int num2 = Int_uniform_ab(0, POPSIZE - 1);
		double f = max(tPopulation[num1].fitness, tPopulation[num2].fitness);//Ҫ��������������нϴ����Ӧ��ֵ
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
 
void elitist()//����ǰ����õģ�����Ӧ��ֵ����
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
 
void evaluate()//�˴���Ӧ��ֵ�����ˮƽ�ߵĴ�С�ĵ�����Ĭ�ϴ����㣬��û���Ų��µ����
{
	//dYTmpÿһ�����һ�Σ����¹�����һ�����ֵ + gap
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

 //���ĵ��ж�ȡ�и���Ϣ�����ɳ�ʼ��Ⱥ
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
	//����sort����Ĭ�����򣬵��ڽṹ����������С�ںţ����������ǽ��򣬼���������������
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

//������ԣ�����n1  n2  n3���������������n1��n2λ�û�����n3ȡ��
void mutate()
{
	const double a = 0.0;
	const double b = 1.0;

	//double w = get_variation_coefficient();
	//cout << "w: " << w << endl;
	//const double k = 0.8;//0 <= k <= 1
	//double dMutateP = k / (5 * (1 + (exp(w))));

	double f_max = 0.0;//��¼��Ⱥ�����Ӧֵ
	double f_sum = 0.0;//��¼��Ⱥ��Ӧ���ܺͣ����Լ�����Ⱥ��Ӧ��ƽ��ֵ
	for(int i = 0; i < POPSIZE; ++i)
	{
		f_max = max(f_max, tPopulation[i].fitness);
		f_sum += tPopulation[i].fitness;
	}
	double f_avg = f_sum / (double)(POPSIZE);
		
	double p_m;//�������
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
	//rand() / double(RAND_MAX)��������0~1֮��ĸ�����
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
		double p = Dou_uniform_ab(a, b);//������0,1������һ��double������
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
	//��ʱ���ʽ��
	len = strftime (time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);
 
	cout << time_buffer << "\n";
 
	return;
# undef TIME_SIZE
}
 
void Xover (int one, int two)
{
	int nPoint = Int_uniform_ab(0, smallRectsNum - 1);

	//�����Ļ����ǲ��Ϸ���
	int * pInt2 = new int [nPoint];//��������󣬴�Ż���2�ظ�������
	int * pInt1 = new int [nPoint];//��������󣬴�Ż���1�ظ�������

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

	//�ѻ���2���ظ����е�����λ�ü�¼���������θ�ֵ������1�ظ����д�
	int * pIntRep2 = new int [cnt];//��Ż���2�ظ������У����ڻ���2�г���˳��
	int * pIntRep1 = new int [cnt];//��Ż���1�ظ������У����ڻ���1�г���˳��

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

		//�޸Ĳ��Ϸ�����1
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

	//�޸Ĳ��Ϸ�����2
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

void Xover2 (int one, int two)//���㽻��
{
	int nPoint1 = Int_uniform_ab(0, smallRectsNum - 1);
	int nPoint2 = Int_uniform_ab(0, smallRectsNum - 1);

	if(nPoint1 > nPoint2)
	{
		int t = nPoint1;
		nPoint1 = nPoint2;
		nPoint2 = t;
	}

	//�����Ļ����ǲ��Ϸ���
	int * pInt2 = new int [nPoint2 - nPoint1];//��������󣬴�Ż���2�ظ�������
	int * pInt1 = new int [nPoint2 - nPoint1];//��������󣬴�Ż���1�ظ�������

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

	//�ѻ���2���ظ����е�����λ�ü�¼���������θ�ֵ������1�ظ����д�
	int * pIntRep2 = new int [cnt];//��Ż���2�ظ������У����ڻ���2�г���˳��
	int * pIntRep1 = new int [cnt];//��Ż���1�ظ������У����ڻ���1�г���˳��

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

		//�޸Ĳ��Ϸ�����1
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

	//�޸Ĳ��Ϸ�����2
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
