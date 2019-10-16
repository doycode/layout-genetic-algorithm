#ifndef __LAYOUT_GA_H__
#define __LAYOUT_GA_H__

typedef struct MYPOINT
{
	double dX;
	double dY;
	MYPOINT(){}

	MYPOINT(double dX_, double dY_)
	{
		dX = dX_;
		dY = dY_;
	}
	//按照y值排序，y相等时按x排
	friend bool operator < (const MYPOINT & tA, const MYPOINT & tB)
	{
		if(tA.dY == tB.dY)
		{
			return tA.dX < tB.dX;
		}
		else
		{
			return tA.dY < tB.dY;
		}
	}
}MYPOINT;

typedef struct MYRECT
{
	double dWidth;
	double dHeight;
	double dS;
	double dRatio;
	int nIndex;
	MYRECT(){}
	MYRECT(double dWidth_, double dHeight_, int nIndex_)
	{
		dWidth = dWidth_;
		dHeight = dHeight_;
		dS = dWidth_ * dHeight_;
		//dRatio = max(dHeight / dWidth, dWidth / dHeight);
		dRatio = (dHeight / dWidth) > (dWidth / dHeight) ? (dHeight / dWidth) : (dWidth / dHeight);
		nIndex = nIndex_;
	}
	MYRECT(double dWidth_, double dHeight_)
	{
		dWidth = dWidth_;
		dHeight = dHeight_;
		dS = dWidth_ * dHeight_;
		//dRatio = max(dHeight / dWidth, dWidth / dHeight);
		dRatio = (dHeight / dWidth) > (dWidth / dHeight) ? (dHeight / dWidth) : (dWidth / dHeight);
	}
	
    friend bool operator < (const MYRECT & tA , const MYRECT & tB)
	{
        double dPriorityA = 0.9 * tA.dS + 0.1 * tA.dRatio;
        double dPriorityB = 0.9 * tB.dS + 0.1 * tB.dRatio;
        return dPriorityA > dPriorityB;//大的靠前
    }
}MYRECT;

typedef struct INDIVIDUAL
{
	int nGene[111];//有正负，负表示矩形旋转，没有0，即编号从1开始，基因编码方式采取矩形编号
    MYPOINT tPosition[111];//每个矩形的左下角位置
	double fitness;
	double rfitness;
	double cfitness;
    INDIVIDUAL(){}
	friend bool operator < (const INDIVIDUAL & tA, const INDIVIDUAL & tB)
	{
		return tA.fitness > tB.fitness;//适应值就是最高水平线的倒数，适应值大的靠前
	}
}INDIVIDUAL;

#if 0
class CLayoutGA
{
public:
	//种群数目，最大迭代次数，单点交叉概率，两点交叉概率，单点旋转变异概率，两点交换变异概率
	CLayoutGA(int, int, double, double, double, double);
	~CLayoutGA();
    void setBigRect(const double dW, const double dH);
    void setGap(const double);
    double getGap();
    void setMargin(const double);
    double getMargin();
    void addRect(const double dW, const double dH, const int nIndex);
    //void setRect(int nPosition, MYRECT tRect);
    MYRECT getRect(int nPosition);
    void dealWithGapMargin();
    void solve();
    void setRectsNum(int);
    int getRectsNum();
    INDIVIDUAL getResult();

public:
	void addRectsToBeCut(int nIndex, MYRECT tRect);
	double getBigRectWidth();
	double getBigRectHeight();

private:
    //最低水平线算法
    double m_dTotalS;
    int m_nSmallRectsNum;//num_of_small_rect
    INDIVIDUAL m_tPopulation[120];
    MYRECT m_tRects[100];//存储要切割的矩形毛坯信息
    double m_dGap;//矩形之间的间距
    double m_dMargin;//与边框的间距
    void initialize();//生成初始种群按照排序规则初始化后存到population里,并更新fitness值
    int* randPerm(int n);//生成N个不同的1-N的随机数
    void fitness(INDIVIDUAL &);
    void crossOne(const INDIVIDUAL tIn1, const INDIVIDUAL tIn2, INDIVIDUAL & tOut1, INDIVIDUAL & tOut2);//单点交叉
    void crossTwo(const INDIVIDUAL tIn1, const INDIVIDUAL tIn2, INDIVIDUAL & tOut1, INDIVIDUAL & tOut2);//双点交叉
    void changePosition(const INDIVIDUAL tIn, INDIVIDUAL & tOut);//位置变异
    void changeRotate(const INDIVIDUAL tIn, INDIVIDUAL & tOut);//旋转变异

private:
	double m_dBigRectWidth;
	double m_dBigRectHeight;


    int m_nMembers;//种群规模
    int m_nMaxIter;//最大迭代次数
    double m_dCrossP1;//单点交叉概率
    double m_dCrossP2;//两点交叉概率
	double m_dMutateP1;//单点旋转变异概率
    double m_dMutateP2;//两点交换变异概率
};
#endif

#endif