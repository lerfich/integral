//
//  main.cpp
//  3rd
//
//  Created by Никита Папенков on 18.05.2020.
//  Copyright © 2020 Никита Папенков. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
using namespace std;

const long double eps = 0.0001;
int rang = 0;
int flag = 0;
int trans = 0;

void printb(vector<double> b){
    for(int i = 0; i < b.size(); ++i)
        cout << b[i] << ' ' << endl;
    
}

void show(vector <vector <double>> A, int n)
{
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cout <<"\t"<< setw(12) << A[i][j] << "\t";
        }
        cout << endl;
    }
}

double NormOfMatrix(vector <vector <double>> A, int n){
    double norma = 0;
    double sum;
    for(int i = 0; i < n; ++i)
    {
       sum = 0;
       for(int j = 0; j < n; ++j)
          sum += abs(A[j][i]);
       if(sum > norma)
           norma = sum;
    }
    return norma;
}

void showvector(vector <double> A, int n)
{
    for(int i = 0; i < n; i++)
    {
        cout << setw(12) << A[i] << "\t";
        cout << endl;
    }
}

vector<vector <double>> MM(vector <vector <double>> A, vector <vector <double>> B, int n)
{
    vector <vector <double>> result(n, vector<double> (n, 0));
    
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            for(int k = 0; k < n; k++)
                result[i][j] += A[i][k] * B[k][j];
            
    /*for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if(abs(result[i][j])<0.00000001) result[i][j] = 0;*/
    
    return result;
}

vector<double> MultMatVec(vector <vector <double>> A, vector<double> x, int n){
    vector<double> result(n, 0);
    double s = 0;
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
                s += A[i][j]*x[j];
    result[i] = s;
    s = 0;
    }
    return result;
}

vector<vector <double>> transpose(vector< vector<double>> A, int n)
{
    double t;
    for(int i = 0; i < n; ++i)
    {
        for(int j = i; j < n; ++j)
        {
            t = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = t;
        }
    }
    return A;
}

vector<double> Solve_Gauss(vector <vector <double>> A, vector <vector <double>> B, vector <double> b, int n, long long int& operat){
    
    for(int i=0; i<n; i++)
    {
        for(int j=i+1;j<n; ++j){
            b[j]-=A[j][i]*b[i];
            operat += 2;
        }
    }
    int i = n;
    while(i--){
        for(int p = n-1; p>i; --p){
            b[i]-=B[i][p]*b[p];
            operat += 2;
        }
        b[i]/=B[i][i];
        ++operat;
    }
    return b;
}

/*std::pair<int, int> maximumMatrixElement(vector<vector<double>> t, int z1, int z2, int f1, int f2){
    double maximum = abs(t[0][0]);
    int a = 0, b = 0;
    for (int i = f1; i < z1; ++i) {
        for (int j = f2; j < z2; ++j) {
            if (abs(t[i][j]) > maximum) {
                maximum = abs(t[i][j]);
                a = i, b = j;
            }
        }
    }
    return std::make_pair(a, b);
}

// LU разложение ведущим элементом
void LUdecomp(vector<vector<double>> A, vector<vector<double>>& L, vector<vector<double>>& U, vector<vector<double>>& Per, vector<vector<double>>& Qer, int n){
    vector<vector<double>> tmp(n, vector<double>(0));
    vector<vector<double>> P(n, vector<double>(0));
    vector<vector<double>> Q(n, vector<double>(0));

    
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++){
            tmp[i].push_back(0);
            P[i].push_back(0);
            Q[i].push_back(0);
        }
    
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if(i == j){P[i][j] = 1; Q[i][j] = 1;}
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            tmp[i][j] = A[i][j];
        
    
    for (int i = 0; i < n - 1; ++i) {
        std::pair<int, int> index = maximumMatrixElement(tmp, n, n, i, i);
        for (int j = 0; j < n; ++j) {
            std::swap(P[0][j], P[index.first][j]);
            std::swap(tmp[0][j], tmp[index.first][j]);
        }
        for (int j = 0; j < n; ++j) {
            std::swap(Q[j][0], Q[j][index.second]);
            std::swap(tmp[j][0], tmp[j][index.second]);
        }
        for (int p = i + 1; p < n; ++p) {
            L[p][i] = tmp[p][i] / tmp[i][i];
            for (int t = i; t < n; ++t) {
                tmp[p][t] = tmp[p][t] - tmp[i][t] * L[p][i];
            }
        }
    }
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            U[i][j] = tmp[i][j];
    
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j){
            Per[i][j] = P[i][j];
            Qer[i][j] = Q[i][j];
        }
}*/

void maxi( vector <vector<double> > A, int currentRow, vector<int>& ret, int dim,  long long int& operat) {
    double maximum = 0;
    for (int i = currentRow; i < dim; ++i) {
        for (int j = currentRow; j <dim; ++j) {
            if (abs(A[i][j]) > maximum) {
                maximum = abs(A[i][j]);
                ret[0] = i;
                ret[1] = j;
                operat += 3;
            }
        }
    }
}


//перестановка строк
void swapRows(vector <vector<double> >& A, int number1, int number2) {
    vector<double> swap;
    swap = A[number1];
    A[number1] = A[number2];
    A[number2] = swap;
}

//перестановка столбцов
void swapColumns(vector <vector<double> >& A, int number1, int number2, int dim,  long long int& operat) {
    double swap;
    for (int i = 0; i < dim; ++i) {
        swap = A[i][number1];
        A[i][number1] = A[i][number2];
        A[i][number2] = swap;
        operat += 3;
    }
}

//перемещение максимального элемента в текущую строку и столбец
void swap_max(vector <vector<double> >& A, vector <vector<double> >& L , vector <vector<double> >& U, vector <vector<double> >& PMatrix, vector <vector<double> >& QMatrix, int currentRow, int n, long long int& operat) {
    vector<int> position_max(2);
    maxi(A, currentRow, position_max, n, operat);
    if (flag == 0) {
        swapRows(A, position_max[0], currentRow);
        swapColumns(A, position_max[1], currentRow, n, operat);
        swapRows(L, position_max[0], currentRow);
        swapColumns(L, position_max[1], currentRow, n, operat);
        swapRows(U, position_max[0], currentRow);
        swapColumns(U, position_max[1], currentRow, n, operat);

        swapRows(PMatrix, position_max[0], currentRow); if (currentRow != position_max[0]) ++trans;
        swapColumns(QMatrix, position_max[1], currentRow, n, operat); if (currentRow != position_max[1]) ++trans;
        operat += 12; //за каждый свап строк +3 операции
        
    }
}


//LU разложение
void LU(vector <vector<double> > A, vector <vector<double> >& L, vector <vector<double> >& U, vector <vector<double> >& P, vector <vector<double> >& Q, int dim, long long int& operat){

    
    for (int i = 0; i < dim; ++i) {
        swap_max(A, L, U, P, Q, i, dim,operat);
        double tmp;
        ++operat;
        L[i][i] = 1;
        for (int j = i + 1; j < dim; ++j) {
            //if(A[i][i] == 0) tmp = 0;
            ++operat;
            tmp = A[j][i] / A[i][i];
            for (int k = i; k < dim; ++k) {
                A[j][k] -= tmp * A[i][k];
                U[j][k] = A[j][k];
                operat+=2;
            }
            A[j][i] = tmp;
            L[j][i] = tmp;
            operat+=2;
        }
    }
    for (int i = 0; i < dim; ++i){
        U[0][i] = A[0][i];
        ++operat;
    }
}


vector <vector <double>> JacobiMatrix(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10){
    vector <vector <double>> J(10, vector<double> (0));
    for(int i = 0;i < 10; ++i)
        for(int j = 0; j < 10; ++j)
            J[i].push_back(0);
    
    J[0][0] = -x2 * sin(x2 * x1);
    J[0][1] = -x1 * sin(x2 * x1);
    J[0][2] = 3.0 * exp(- (3.0 * x3));
    J[0][3] = x5 * x5;
    J[0][4] = 2.0 * x4 * x5;
    J[0][5] = -1.0;
    J[0][6] = 0.0;
    J[0][7] = -2.0 * cosh(2.0 * x8) * x9;
    J[0][8] = -sinh(2.0 * x8);
    J[0][9] = 2.0;
    J[1][0] = x2 * cos(x2 * x1);
    J[1][1] = x1 * cos(x2 * x1);
    J[1][2] = x9 * x7;
    J[1][3] = 0.0;
    J[1][4] = 6.0 * x5;
    J[1][5] = -exp(-x10 + x6) -  x8 - 1.0;
    J[1][6] =  x3 * x9;
    J[1][7] = -x6;
    J[1][8] =  x3 * x7;
    J[1][9] = exp(-x10 + x6);
    J[2][0] = 1;
    J[2][1] = -1;
    J[2][2] = 1;
    J[2][3] = -1;
    J[2][4] = 1;
    J[2][5] = -1;
    J[2][6] = 1;
    J[2][7] = -1;
    J[2][8] = 1;
    J[2][9] = -1;
    J[3][0] = - x5 * pow( x3 + x1, -2.0);
    J[3][1] = -2.0 * x2 * cos(x2 * x2);
    J[3][2] = - x5 * pow( x3 + x1, -2.0);
    J[3][3] = -2.0 * sin(-x9 +  x4);
    J[3][4] = 1.0 / ( x3 + x1);
    J[3][5] = 0;
    J[3][6] = -2.0 * cos(x7 * x10) * x10 * sin(x7 * x10);
    J[3][7] = -1;
    J[3][8] = 2.0 * sin(-x9 +  x4);
    J[3][9] = -2.0 * cos(x7 * x10) * x7 * sin(x7 * x10);
    J[4][0] = 2 * x8;
    J[4][1] = -2.0 * sin(x2);
    J[4][2] = 2 * x8;
    J[4][3] = pow(-x9 +  x4, -2.0);
    J[4][4] = cos( x5);
    J[4][5] = x7 * exp(-x7 * (-x10 + x6));
    J[4][6] = -(x10 - x6) * exp(-x7 * (-x10 + x6));
    J[4][7] =  (2 * x3) + 2.0 * x1;
    J[4][8] = -pow(-x9 +  x4, -2.0);
    J[4][9] = -x7 * exp(-x7 * (-x10 + x6));
    J[5][0] = exp(x1 -  x4 - x9);
    J[5][1] = -3.0 / 2.0 * x10 * sin(3.0 * x10 * x2);
    J[5][2] = -x6;
    J[5][3] = -exp(x1 -  x4 - x9);
    J[5][4] = 2 * x5 / x8;
    J[5][5] = -x3;
    J[5][6] = 0;
    J[5][7] = -x5 * x5 * pow( x8,  (-2));
    J[5][8] = -exp(x1 -  x4 - x9);
    J[5][9] = -3.0 / 2.0 * x2 * sin(3.0 * x10 * x2);
    J[6][0] = cos( x4);
    J[6][1] = 3.0 * x2 * x2 * x7;
    J[6][2] = 1;
    J[6][3] = -(x1 - x6) * sin( x4);
    J[6][4] = x10 * pow( x5,  (-2)) * cos(x10 /  x5 +  x8);
    J[6][5] = -cos( x4);
    J[6][6] = pow(x2, 3.0);
    J[6][7] = -cos(x10 /  x5 +  x8);
    J[6][8] = 0;
    J[6][9] = -1.0 /  x5 * cos(x10 /  x5 +  x8);
    J[7][0] = 2.0 *  x5 * (x1 - 2.0 * x6);
    J[7][1] = -x7 * exp(x2 * x7 + x10);
    J[7][2] = -2.0 * cos(-x9 +  x3);
    J[7][3] = 1.5;
    J[7][4] = pow(x1 - 2.0 * x6, 2.0);
    J[7][5] = -4.0 *  x5 * (x1 - 2.0 * x6);
    J[7][6] = -x2 * exp(x2 * x7 + x10);
    J[7][7] = 0;
    J[7][8] = 2.0 * cos(-x9 +  x3);
    J[7][9] = -exp(x2 * x7 + x10);
    J[8][0] = -3;
    J[8][1] = -2.0 *  x8 * x10 * x7;
    J[8][2] = 0;
    J[8][3] = exp( (x5 + x4));
    J[8][4] = exp( (x5 + x4));
    J[8][5] = -7.0 * pow(x6, -2.0);
    J[8][6] = -2.0 * x2 *  x8 * x10;
    J[8][7] = -2.0 * x2 * x10 * x7;
    J[8][8] = 3;
    J[8][9] = -2.0 * x2 *  x8 * x7;
    J[9][0] = x10;
    J[9][1] = x9;
    J[9][2] = -x8;
    J[9][3] = cos( x4 +  x5 + x6) * x7;
    J[9][4] = cos( x4 +  x5 + x6) * x7;
    J[9][5] = cos( x4 +  x5 + x6) * x7;
    J[9][6] = sin( x4 +  x5 + x6);
    J[9][7] = -x3;
    J[9][8] = x2;
    J[9][9] = x1;
    return J;
}

vector<double> SU(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10){
    vector <double> F(10, 0);
    for(int i = 0; i < 10; ++i)
        F.push_back(0.0);
    F[0] = cos(x2 * x1) - exp(- (3.0 * x3)) + x4 * x5 * x5 - x6 - sinh( (2.0 * x8)) * x9 +  (2.0 * x10) + 2.000433974165385440;
    F[1] = sin(x2 * x1) +  x3 * x9 * x7 - exp(- x10 + x6) + 3.0 * x5 * x5 - x6 *  (x8 + 1.0) + 10.886272036407019994;
    F[2] = x1 - x2 +  x3 - x4 + x5 - x6 + x7 -  x8 + x9 -  x10 - 3.1361904761904761904;
    F[3] = 2.0 * cos(-x9 + x4) + x5 / ( x3 + x1) - sin(x2 * x2) + pow(cos(x7 *  x10), 2.0) -  x8 - 0.1707472705022304757;
    F[4] = sin(x5) + 2.0 *  x8 * ( x3 + x1) - exp(-x7 * (- x10 + x6)) + 2.0 * cos(x2) - 1.0 / (-x9 + x4) - 0.3685896273101277862;
    F[5] = exp(x1 - x4 - x9) + x5 * x5 /  x8 + cos(3.0 *  x10 * x2) / 2.0 - x6 *  x3 + 2.0491086016771875115;
    F[6] = pow(x2, 3.0) * x7 - sin( x10 / x5 +  x8) + (x1 - x6) * cos(x4) +  x3 - 0.7380430076202798014;
    F[7] = x5 * pow(x1 - 2.0 * x6, 2.0) - 2.0 * sin(-x9 +  x3) + 1.5 * x4 - exp(x2 * x7 +  x10) + 3.5668321989693809040;
    F[8] = 7.0 / x6 + exp(x5 + x4) - 2.0 * x2 *  x8 *  x10 * x7 + 3.0 * x9 - 3.0 * x1 - 8.4394734508383257499;
    F[9] =  x10 * x1 + x9 * x2 -  (x8 * x3) + sin(x4 + x5 + x6) * x7 - 0.78238095238095238096;
    F[0] *= -1; F[1] *= -1; F[2] *= -1; F[3] *= -1; F[4] *= -1; F[5] *= -1; F[6] *= -1; F[7] *= -1; F[8] *= -1; F[9] *= -1;
    return F;
}

vector<double> SP(vector<double> x0){
    x0[0] = 0.5;
    x0[1] = 0.5;
    x0[2] = 1.5;
    x0[3] = -1.0;
    x0[4] = -0.5;
    x0[5] = 1.5;
    x0[6] = 0.5;
    x0[7] = -0.5;
    x0[8] = 1.5;
    x0[9] = -1.5;
    return x0;
}

bool OK(vector<double> dx, double eps){
    bool f = true;
    for(int i = 0; i < 10; ++i)
        if(dx[i] >= eps)
            f = false;
    return f;
}

void Newton(vector<double> x, vector<double> x0){
    
    clock_t start = clock();
    long long int iterat = 0;
    long long int operat = 0;
    vector<double> dx(10, 0.0);
    vector<double> NewF(10, 0.0);
    
    vector <vector <double>> L(10, vector<double> (0));
    vector <vector <double>> U(10, vector<double> (0));
    vector <vector <double>> P(10, vector<double> (0));
    vector <vector <double>> Q(10, vector<double> (0));
    vector <vector <double>> NewJac(10, vector<double> (0));
    
    for(int i = 0; i < 10; ++i)
        dx.push_back(0.0);
    
    for(int i = 0; i < 10; ++i){
        for(int j = 0; j < 10; ++j){
            NewJac[i].push_back(0);
            L[i].push_back(0);
            U[i].push_back(0);
            P[i].push_back(0);
            Q[i].push_back(0);
            if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
        }
        NewF.push_back(0);
    }
    
    for(int i = 0; i < 10; ++i)
        x[i] = x0[i];
    
    do{
        NewJac = JacobiMatrix(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        LU(NewJac, L, U, P, Q, 10, operat);
        NewF = SU(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        dx = MultMatVec(Q,Solve_Gauss(L, U, MultMatVec(P, NewF, 10), 10, operat),10);
        for(int i = 0; i < 10; ++i)
            x[i] += dx[i];
        ++iterat;
        
        for(int i = 0; i < 10; ++i)
            for(int j = 0; j < 10; ++j){
                L[i][j] = 0; U[i][j] = 0; P[i][j] = 0; Q[i][j] = 0;
                if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
            }
        
    }while(!OK(dx, eps));
    double sec = ((double)(clock() - start)) / CLOCKS_PER_SEC ;
    cout << "    X with the Newton method:    " << endl;
    for(int i = 0; i < 10; ++i)
        cout << setw(15) << x[i] << endl;
    cout <<  " Iterations : " << iterat << endl;
    cout <<  " Arithmetic operations : " << operat << endl;
    cout <<  " Time : " << sec << endl;
}

void ModifiedNewton(vector<double> x, vector<double> x0){
    clock_t start = clock();
    long long int iterat = 0;
    long long int operat = 0;
    vector<double> dx(10, 0.0);
    vector<double> NewF(10, 0.0);
    
    vector <vector <double>> L(10, vector<double> (0));
    vector <vector <double>> U(10, vector<double> (0));
    vector <vector <double>> P(10, vector<double> (0));
    vector <vector <double>> Q(10, vector<double> (0));
    vector <vector <double>> NewJac(10, vector<double> (0));
    
    for(int i = 0; i < 10; ++i)
        dx.push_back(0.0);
    
    for(int i = 0; i < 10; ++i){
        for(int j = 0; j < 10; ++j){
            NewJac[i].push_back(0);
            L[i].push_back(0);
            U[i].push_back(0);
            P[i].push_back(0);
            Q[i].push_back(0);
            if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
        }
        NewF.push_back(0);
    }
    
    for(int i = 0; i < 10; ++i)
        x[i] = x0[i];
    
    NewJac = JacobiMatrix(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
    LU(NewJac, L, U, P, Q, 10, operat);
    operat = 0;
    
    do{
        NewF = SU(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        dx = MultMatVec(Q,Solve_Gauss(L, U, MultMatVec(P, NewF, 10), 10, operat),10);
        for(int i = 0; i < 10; ++i)
            x[i] += dx[i];
        ++iterat;
    }while(!OK(dx, eps));
    double sec = ((double)(clock() - start)) / CLOCKS_PER_SEC ;
    cout << "    X with the Modified Newton method:    " << endl;
    for(int i = 0; i < 10; ++i)
        cout << setw(15) << x[i] << endl;
    cout <<  " Iterations : " << iterat << endl;
    cout <<  " Arithmetic operations : " << operat << endl;
    cout <<  " Time : " << sec << endl;
}

void SwitchNM(vector<double> x, vector<double> x0, int k){
    clock_t start = clock();
    long long int t = 0;
    long long int iterat = 0;
    long long int operat = 0;
    vector<double> dx(10, 0.0);
    vector<double> NewF(10, 0.0);
    
    vector <vector <double>> L(10, vector<double> (0));
    vector <vector <double>> U(10, vector<double> (0));
    vector <vector <double>> P(10, vector<double> (0));
    vector <vector <double>> Q(10, vector<double> (0));
    vector <vector <double>> NewJac(10, vector<double> (0));
        
    for(int i = 0; i < 10; ++i)
        dx.push_back(0.0);
    for(int i = 0; i < 10; ++i){
        for(int j = 0; j < 10; ++j){
            NewJac[i].push_back(0);
            L[i].push_back(0);
            U[i].push_back(0);
            P[i].push_back(0);
            Q[i].push_back(0);
            if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
        }
        NewF.push_back(0);
    }
        
    for(int i = 0; i < 10; ++i)
        x[i] = x0[i];
    while(k>0){
        NewJac = JacobiMatrix(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        LU(NewJac, L, U, P, Q, 10, operat);
        NewF = SU(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        dx = MultMatVec(Q,Solve_Gauss(L, U, MultMatVec(P, NewF, 10), 10, operat),10);
        for(int i = 0; i < 10; ++i)
            x[i] += dx[i];
        ++iterat;
        for(int i = 0; i < 10; ++i)
            for(int j = 0; j < 10; ++j){
                L[i][j] = 0; U[i][j] = 0; P[i][j] = 0; Q[i][j] = 0;
                if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
            }
        --k;
    }
    
    t = operat;
    NewJac = JacobiMatrix(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
    LU(NewJac, L, U, P, Q, 10, operat);
    operat = t;
    
    do{
        NewF = SU(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        dx = MultMatVec(Q,Solve_Gauss(L, U, MultMatVec(P, NewF, 10), 10, operat),10);
        for(int i = 0; i < 10; ++i)
            x[i] += dx[i];
        ++iterat;
    }while(!OK(dx, eps));
    double sec = ((double)(clock() - start)) / CLOCKS_PER_SEC ;
    cout << "    X counted with switching between Newton and Modified Newton methods:    " << endl;
    for(int i = 0; i < 10; ++i)
        cout << setw(15) << x[i] << endl;
    cout <<  " Iterations : " << iterat << endl;
    cout <<  " Arithmetic operations : " << operat << endl;
    cout <<  " Time : " << sec << endl;
}

void HybridNewton(vector<double> x, vector<double> x0, int k){
    clock_t start = clock();
    long long int t = 0;
    long long int iterat = 0;
    long long int operat = 0;
    vector<double> dx(10, 0.0);
    vector<double> NewF(10, 0.0);
    
    vector <vector <double>> L(10, vector<double> (0));
    vector <vector <double>> U(10, vector<double> (0));
    vector <vector <double>> P(10, vector<double> (0));
    vector <vector <double>> Q(10, vector<double> (0));
    vector <vector <double>> NewJac(10, vector<double> (0));
    
    for(int i = 0; i < 10; ++i)
        dx.push_back(0.0);
    
    for(int i = 0; i < 10; ++i){
        for(int j = 0; j < 10; ++j){
            NewJac[i].push_back(0);
            L[i].push_back(0);
            U[i].push_back(0);
            P[i].push_back(0);
            Q[i].push_back(0);
            if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
        }
        NewF.push_back(0);
    }
    
    for(int i = 0; i < 10; ++i)
        x[i] = x0[i];
    
    do{
        if(!(iterat%k)){
            for(int i = 0; i < 10; ++i)
                for(int j = 0; j < 10; ++j){
                    L[i][j] = 0; U[i][j] = 0; P[i][j] = 0; Q[i][j] = 0;
                    if(i == j) {P[i][j] = 1; Q[i][j] = 1;}
                }
            t = operat;
            NewJac = JacobiMatrix(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
            LU(NewJac, L, U, P, Q, 10, operat);
            operat = t;
        }
        NewF = SU(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
        dx = MultMatVec(Q,Solve_Gauss(L, U, MultMatVec(P, NewF, 10), 10, operat),10);
        for(int i = 0; i < 10; ++i)
            x[i] += dx[i];
        ++iterat;
    }while(!OK(dx, eps));
    double sec = ((double)(clock() - start)) / CLOCKS_PER_SEC ;
    cout << "    X counted with Hybrid Newton method:    " << endl;
    for(int i = 0; i < 10; ++i)
        cout << setw(15) << x[i] << endl;
    cout <<  " Iterations : " << iterat << endl;
    cout <<  " Arithmetic operations : " << operat << endl;
    cout <<  " Time : " << sec << endl;
}


int main(){
    int n = 10;
    int k;
    cout << "Как часто пересчитывать матрицу Якоби для гибридного метода Ньютона?" << endl;
    cin >> k;
    vector<double> F(n, 0), x0(n,0), x(n,0);
    vector <vector <double>> A(n, vector<double> (0)), J(n, vector<double> (0));/* L(n, vector<double> (0)), U(n, vector<double> (0)), P(n, vector<double> (0)), Q(n, vector<double> (0));*/
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i].push_back(0);
            J[i].push_back(0);
        }
        F.push_back(0);
        x0.push_back(0);
        x.push_back(0);
    }
    //JacobiMatrix(J, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
    //F = SU(F, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9]);
    x0=SP(x0);
    Newton(x, x0);
    cout << endl;
    ModifiedNewton(x, x0);
    cout << endl;
    SwitchNM(x, x0, 4);
    cout << endl;
    HybridNewton(x, x0, k);
    cout << endl;
    
    x0[4] = -0.2;
    
    cout << " After changing x5 from -0.5 to -0.2 " << endl;
    cout << endl;
    
    Newton(x, x0);
    cout << endl;
    /*ModifiedNewton(x, x0);
    cout << endl;*/
    SwitchNM(x, x0, 8); // при k < 7, либо ломается, либо большая ошибка, при k = 7, близкий ответ, при k > 7 точный ответ
    cout << endl;
    HybridNewton(x, x0, k);
    return 0;
}


// 5               k=3
// 7685
// 0.001535


// 6               k=2
// 5644
// 0.002354


// 5               k=4
// 9931
// 0.001218


// 6               k=5
// 12367
// 0.001614
