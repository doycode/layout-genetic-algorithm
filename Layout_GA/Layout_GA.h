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
	//����yֵ����y���ʱ��x��
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
        return dPriorityA > dPriorityB;//��Ŀ�ǰ
    }
}MYRECT;

typedef struct INDIVIDUAL
{
	int nGene[111];//������������ʾ������ת��û��0������Ŵ�1��ʼ��������뷽ʽ��ȡ���α��
    MYPOINT tPosition[111];//ÿ�����ε����½�λ��
	double fitness;
	double rfitness;
	double cfitness;
    INDIVIDUAL(){}
	friend bool operator < (const INDIVIDUAL & tA, const INDIVIDUAL & tB)
	{
		return tA.fitness > tB.fitness;//��Ӧֵ�������ˮƽ�ߵĵ�������Ӧֵ��Ŀ�ǰ
	}
}INDIVIDUAL;

#if 0
class CLayoutGA
{
public:
	//��Ⱥ��Ŀ�����������������㽻����ʣ����㽻����ʣ�������ת������ʣ����㽻���������
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
    //���ˮƽ���㷨
    double m_dTotalS;
    int m_nSmallRectsNum;//num_of_small_rect
    INDIVIDUAL m_tPopulation[120];
    MYRECT m_tRects[100];//�洢Ҫ�и�ľ���ë����Ϣ
    double m_dGap;//����֮��ļ��
    double m_dMargin;//��߿�ļ��
    void initialize();//���ɳ�ʼ��Ⱥ������������ʼ����浽population��,������fitnessֵ
    int* randPerm(int n);//����N����ͬ��1-N�������
    void fitness(INDIVIDUAL &);
    void crossOne(const INDIVIDUAL tIn1, const INDIVIDUAL tIn2, INDIVIDUAL & tOut1, INDIVIDUAL & tOut2);//���㽻��
    void crossTwo(const INDIVIDUAL tIn1, const INDIVIDUAL tIn2, INDIVIDUAL & tOut1, INDIVIDUAL & tOut2);//˫�㽻��
    void changePosition(const INDIVIDUAL tIn, INDIVIDUAL & tOut);//λ�ñ���
    void changeRotate(const INDIVIDUAL tIn, INDIVIDUAL & tOut);//��ת����

private:
	double m_dBigRectWidth;
	double m_dBigRectHeight;


    int m_nMembers;//��Ⱥ��ģ
    int m_nMaxIter;//����������
    double m_dCrossP1;//���㽻�����
    double m_dCrossP2;//���㽻�����
	double m_dMutateP1;//������ת�������
    double m_dMutateP2;//���㽻���������
};
#endif

#endif