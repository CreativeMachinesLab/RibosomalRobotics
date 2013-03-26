/* ---------------------------------------------------
   FILE:     bodyPlan.h
	AUTHOR:   Josh Bongard
	DATE:     October 2, 2000
	FUNCTION: This class contains all information for
				 a single physical segment of an organism
				 in the MathEngine environment.
 -------------------------------------------------- */
#include "iostream"
#include "fstream"

using namespace std;

#ifndef _MATRIX_H
#define _MATRIX_H

class MATRIX {

public:
	int length;
	int width;
	double *vals;

public:
	MATRIX(int ln, int wd);
	MATRIX(int ln, int wd, double val);
	MATRIX(MATRIX *m);
	MATRIX(ifstream *inFile);
	MATRIX(int ln, int wd, ifstream *inFile, int dummyVariable);
	MATRIX(char *fileName);
	~MATRIX(void);
	double   AbsoluteDifference(MATRIX *m);
	double   AbsoluteDifference(MATRIX *m, int startWd, int endWd);
	double   AbsoluteDifference(MATRIX *m, int wd);
	double   AbsoluteDifference(MATRIX *m, int wd, double normFactor);
	double   AbsoluteDifferenceOfTails(MATRIX *m);
	void	 Add(MATRIX *other);
	void	 Add(int i, int j, double val);
	void     AddToColumn(int j, double val);
	void	 AppendToFile(char *fileName);
	void	 AverageMatrices(MATRIX **matrices, int numMatrices);
	int      BinaryToDecimal(int i, int j, int ln);
	int		 ComputeFSM(MATRIX *str, MATRIX *states, int i, int strIndex, int strLength, int q0, MATRIX *F);
	void     CopyRow(int myRow, MATRIX *other, int hisRow);
	MATRIX  *CountValues(int i1, int i2, int j1, int j2, int maxVal);
	void	 CreateChain(int depth);
	void	 CreateIdentity(void);
	void	 CreateParity(void);
	void     CreatePermutation(int min, int max);
	void	 DecreaseFullColumns(void);
	double   DistBetRows(int i1, int i2);
	void     Div(double val);
	void     Div(int i, int j, double val);
	double   EqualColumnVals(int col1, int col2);
	void	 FillColumn(int myCol, int hisCol, MATRIX *m);
	void     FillRow(int myRow, int hisRow, MATRIX *m);
	void     FindMinMaxInColumn(int j,int i1, int i2,double *min, double *max);
	void     Flip(int i, int j);
	void	 FlipRandomBit(void);
	void	 FlipRandomBitInRow(int i);
	double   Get(int i, int j);
	MATRIX  *GetColumn(int j);
	double   GetDistanceReading(void);
	double   GetForwardBackTiltReading(void);
	double   GetLeastReadingFromPressureSensors(void);
	double   GetLeftRightTiltReading(void);
	double   GetMaxReadingFromPressureSensors(void);
	void     GetMaxValsInColumn(int j, MATRIX *maxVals);
	void     GetMaxValsInColumn(int j, int iFirst, int iLast, MATRIX *maxVals);
	MATRIX  *GetRow(int i);
	double   GetTotalReadingsFromPressureSensors(void);
	double   GetTotalSpeed(void);
	int		 In(double val);
	void     IncreaseEmptyColumns(void);
	void	 InitColumns(int colSum);
	void	 Load(ifstream *inFile);
	double   MaxDiff(MATRIX *otherM);
	double   MaxDistBetRows(MATRIX *otherM);
	int	 MaxIndexInRow(int i);
	double   MaxValInColumn(int j);
	double   MaxValInColumnContingentOn(int j, int k);
	double   Mean(void);
	double   Mean(int i1, int i2, int j1, int j2);
	double   MeanDist(MATRIX *otherM);
	double   MeanDist(MATRIX *otherM, int i1, int i2, int j1, int j2);
	double   MeanDistBetRows(void);
	double	 MeanDistBetRows(MATRIX *otherM);
	double   MeanUpperTriangle(void);
	double	 Min(void);
	double   MinDistBetRows(void);
	double   MinValInColumn(int j);
	int		 MostSimilarRow(MATRIX *r, int i1, int i2);
	double	 MSE(MATRIX *other);
	void	 Mult(double val);
	void     Mult(int i, int j, double val);
	void     MultColumn(int j, double val);
	void     Mutate(void);
	void	 Normalize(void);
	void	 Nudge(double maxVal);
	void     Permute(void);
	void	 Perturb(double maxVal);
	void	 Print(void);
	void	 Print(int i1, int i2, int j1, int j2);
	void	 PrintColumn(int j);
	void     Randomize(int maxVal);
	void     RandomizeColumn(int j, int maxVal);
	void     RandomizeRow(int i, int maxVal);
	void	 Reduce(int amt);
	void	 Replace(MATRIX *m);
	void	 ReZero(void);
	double	 RollingMean(MATRIX *other, int h, int w);
	void	 Round(void);
	int      RowContainsEqualValues(int i);
	int		 RowDifference(int myRow, int hisRow, MATRIX *m);
	void	 SelectUniquelyFrom(int maxVal);
	void	 Set(int i, int j, double val);
	void	 SetAllTo(double val);
	void     SetRow(int i, int val);
	void     Sqrt(void);
	void     Sqrt(int i, int j);
	double   Sum(void);
	double	 Sum(int i1, int i2, int j1, int j2);
	int		 SumOfColumn(int j);
	int		 SumOfColumn(int i1, int i2, int j);
	int		 SumOfIndices(MATRIX *indices,int i, int j1, int j2);
	double	 SumOfRow(int i);
	double   SumOfRow(int i, int j1, int j2);
	void	 SwapRows(MATRIX *other);
	void     SwapValues(MATRIX *other);
	void	 Write(ofstream *outFile);
	void	 WriteAndRename(char *fileName);
	void	 ZeroColumn(int j);

private:
	int		BinaryToDecimal(int j);
	int	    Contains(int val);
	MATRIX *DecimalToBinary(int val, int maxValue);
	int		FindFirstValue(int j);
	double  Mean(int i, int startJ, int endJ);
	void	InitColumn(int j, int colSum);
	void	InitRandomly(void);
	void    ReadFromFile(ifstream *inFile);
	double  RollingMean(MATRIX *other, int j, int h, int w);
};

#endif
