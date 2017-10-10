//
//  main.cpp
//  Test Point Cloud Project
//
//  Created by 李喆 on 2017/9/17.
//  Copyright © 2017年 李喆. All rights reserved.
//

//#include <iostream>
//#include <string>
//using namespace std;
//
//int main(int argc, const char * argv[]) {
//    // insert code here...
//    std::string s;
//    std::cin >> s;
//    std::cout << "Hello, " << s << "!\n";
//    return 0;
//}

//------------------------------------------------------------------------------
//    Copyright (c) 2012 eryar All Rights Reserved.
//
//        File    : Main.cpp
//        Author  : eryar@163.com
//        Date    : 2012-9-9 17:11
//        Version : 0.1v
//
//    Description : lpsolve test program.
//
//==============================================================================
#include <iostream>
#include "stdio.h"
#include <fstream>
#include <vector>
#include <math.h>
using namespace std;

#include "lp_lib.h"

//#pragma comment(lib, "lpsolve55.lib")

int demo(void);
void readPtsFile(string filename);
void readPtsFileCVersion(char* filename);
int countLines(string filename);
int ** creatPointCloudArray(int rows);
vector<vector<double>> creatPointCloudArrayFromFile(string filename);
void printVector(vector<vector<double>> pointclouds);
REAL calculateRotation(int i, int j, int k, int l, const vector<vector<double>> &pointclouds1, const vector<vector<double>> &pointclouds2);
REAL mylpsolve(int m ,int n, vector<vector<double>> &pointclouds1, vector<vector<double>> &pointclouds2, double threshold);

int main(int argc, char* argv[])
{
    vector<vector<double>> pointclouds1 = creatPointCloudArrayFromFile("/Users/lizhe/Desktop/pointclouddataset/test1.pts");
    vector<vector<double>> pointclouds2 = creatPointCloudArrayFromFile("/Users/lizhe/Desktop/pointclouddataset/test2.pts");
    //vector<vector<double>> pointclouds3 = creatPointCloudArrayFromFile("/Users/lizhe/Desktop/pointclouddataset/000020.pts");
    printVector(pointclouds1);
    cout << "========" << endl;
    //printVector(pointclouds2);
    cout << "========" << endl;
    //printVector(pointclouds3);
    //cout << "========" << endl;
    
    time_t start,stop;
    start = time(NULL);
    REAL theta = mylpsolve(pointclouds1.size(), pointclouds2.size(), pointclouds1, pointclouds2, 1);
    stop = time(NULL);
    cout << "theta: " << theta << " ------- time usage: "<< stop-start <<endl;
   
// ========================================
    
//    return demo();

// ========================================
    
//    // the following is some test for point cloud file reading
//    readPtsFile("/Users/lizhe/Desktop/lalala.txt"); //OK
//    readPtsFile("test1.pts"); // need to put the file in the directory of "Products", not the same directory of main
//    char filename[] = "lalala.txt";
//    readPtsFileCVersion(filename);
    
// ========================================
    
//    // the following is some test for point cloud vector storage
//    vector<vector<double>> pointclouds = creatPointCloudArrayFromFile("/Users/lizhe/Desktop/pointclouddataset/test2.pts");
//    vector<vector<double>> pointclouds = creatPointCloudArrayFromFile("/Users/lizhe/Desktop/pointclouddataset/000020.pts");
//    printVector(pointclouds);
   
// ========================================
    
//    // the following is some test for point cloud array storage
//    int rows = countLines("/Users/lizhe/Desktop/lalala.txt");
//
//    int ** pointcloud = creatPointCloudArray(rows);
//    for (int i = 0; i < rows; i++){
//        for (int j = 0; j < 3; j++){
//            cout << pointcloud[i][j] << "  ";
//        }
//        cout << endl;
//    }

// ========================================
    
//    // dynamic create 2D array
//    int **pointcloud = new int*[rows];
//    for (int i = 0; i < rows; i++){
//        pointcloud[i] = new int[3];
//    }
//
//
//    // don't forget to free the memory after usage
//    for (int i = 0; i < rows; i++){
//        delete []pointcloud[i];
//        pointcloud[i] = NULL;
//    }
//    delete []pointcloud;
//    pointcloud = NULL;
    
}

REAL mylpsolve(int m ,int n, vector<vector<double>> &pointclouds1, vector<vector<double>> &pointclouds2, double threshold){
   
    lprec* lp;
    int ret = 0;
    REAL thetaResult = -1;
    
    // the m*n+1 variable is theta
    lp = make_lp(0, m*n+1);
    
    if (NULL == lp){
        cout << "unable to create new LP nodel!\n" << endl;
        ret = 1;
    }

    // sets all variable(except theta) to binary
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            if (!set_binary(lp, i*n+j+1, TRUE)){
                ret = 2;
                cout << "set variable to binary fail" << endl;
            }
        }
    }
    
    // the official site says that declare the objective function than constrain is much better
    if (0 == ret) {
        // Row mode should be truned off again when done building the model.
        set_add_rowmode(lp, FALSE);
        
        int colno[1] = {m*n+1};
        REAL row[1] = {1.0};
        
        // Set the objective in lpsolve.
        if (!set_obj_fnex(lp, 1, row, colno))
        {
            ret = 4;
            cout << "set objective fail" << endl;
        }
    }
    
    // Makes building the model faster if it is done rows by row.
    set_add_rowmode(lp, TRUE);
    
    if (0 == ret) {
        
        // add constraint for Sum(Sik) >= 1, for i belong to m, totally m equations
        // constructs the row like : S11 + S12 + S13  ...  Sin >= 1
        REAL sparserowX[n]; // the coefficient, e.g. 1
        int colnoX[n];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                colnoX[j] = i*n+j+1;
                sparserowX[j] = 1.0;
               //cout << colnoX[j] << " " << sparserowX[j] << endl;
            }
            if(!add_constraintex(lp, n, sparserowX, colnoX, GE, 1.0)){
                ret = 3;
                cout << "add constaint fail in step 1" << endl;
            }
        }
    }
    
    if (0 == ret) {
        // add constraint for Sum(Sik) >= 1, for i belong to n, totally n equations
        // constructs the row like : S11 + S21 + S31  ... Sm1 >= 1
        REAL sparserowY[m];
        int colnoY[m];
        for (int j = 0; j < n; j++){
            for (int i = 0; i < m; i++){
                colnoY[i] = i*n+j+1;
                sparserowY[i] = 1.0;
               // cout << colnoY[i] << " " << sparserowY[i] << endl;
            }
            if (!add_constraintex(lp, m, sparserowY, colnoY, GE, 1.0)) {
                ret = 3;
                cout << "add constaint fail in step 2" << endl;
            }
        }
    }
    
    if (0 == ret) {
        // the colnoXY[0] colnoXY[1] and sparserowXY[2] will be change
        REAL sparserowXY[3] = {1.0,1.0,1.0};
        int colnoXY[3] = {1,1,m*n+1};
        // constructs the row like: Sik + Sjl - theta/distortion <= 1
        for (int i1 = 0; i1 < m; i1++){
            for (int j1 = 0; j1 < n; j1++){
                for (int i2 = 0; i2 < m; i2 ++){
                    for (int j2 = 0; j2 < n; j2++){
                        colnoXY[0] = i1*n + j1 + 1;
                        colnoXY[1] = i2*n + j2 + 1;
                        sparserowXY[2] = -1 / calculateRotation(i1, i2, j1, j2, pointclouds1, pointclouds2);
                        if (!add_constraintex(lp, 3, sparserowXY, colnoXY, LE, 1.0)) {
                            ret = 3;
                            cout << "add constaint fail in step 3" << endl;
                        }
                       // cout << colnoXY[0] << " " << colnoXY[1] << " " << sparserowXY[2] << endl;
                    }
                }
            }
        }
    }
    
    
    // A solution is calculated, now lets get some results.
    if (0 == ret)
    {
        ret = solve(lp);
        
        cout << "ret: " << ret << endl;
        
        // Objective value.
        thetaResult = get_objective(lp);
        
        //cout<<"Objective value: "<< thetaResult <<endl;
        
        REAL * variables = new REAL[m*n+1];
        get_variables(lp, variables);
        for (int i = 0; i < m*n+1; i++)
        {
            cout << variables[i] <<endl;
        }
    }
    
    // Clean up such that all used memory by lpsolve is freed.
    if (lp != NULL)
    {
        delete_lp(lp);
    }
    
    return thetaResult;
}

REAL calculateRotation(int i, int j, int k, int l, const vector<vector<double>> &pointclouds1, const vector<vector<double>> &pointclouds2){
    
    double diffX1 = fabs(pointclouds1[i][0] - pointclouds1[j][0]);
    double diffY1 = fabs(pointclouds1[i][1] - pointclouds1[j][1]);
    double diffZ1 = fabs(pointclouds1[i][2] - pointclouds1[j][2]);
    double diffX2 = fabs(pointclouds2[k][0] - pointclouds2[l][0]);
    double diffY2 = fabs(pointclouds2[k][1] - pointclouds2[l][1]);
    double diffZ2 = fabs(pointclouds2[k][2] - pointclouds2[l][2]);
    
    double distance1 = sqrt(diffX1 * diffX1 + diffY1 * diffY1 + diffZ1 * diffZ1);
    double distance2 = sqrt(diffX2 * diffX2 + diffY2 * diffY2 + diffZ2 * diffZ2);
    
    REAL result = distance1 >= distance2 ? distance1 : distance2;
    return result == 0.0 ? result : 0.000001 ;
}

int demo( void )
{
    lprec*  lp;
    int Ncol    = 0;
    int *colno  = NULL;
    int j       = 0;
    int ret     = 0;
    REAL* row   = NULL;
    
    // We will build the model row by row,
    // So we start with creating a model with 0 rows and 2 columns.
    
    // There are two bariables in the model.
    Ncol    = 2;
    
    lp  = make_lp(0, Ncol);
    
    if (lp == NULL)
    {
        cout<<"Unable to create new LP model!\n"<<endl;
        ret = 1;
    }
    
    if (ret == 0)
    {
        // Let us name our bariables. Not required, but can be useful for debugging.
        
        char xcord[] = "x";
        char ycord[] = "y";
        
        set_col_name(lp, 1, xcord);
        set_col_name(lp, 2, ycord);
        
        // Create space large enough for one row.
        colno   = (int *)malloc(Ncol * sizeof(*colno));
        row     = (REAL*)malloc(Ncol * sizeof(*row));
        
        // malloc memory failed.
        if ((colno == NULL) || row == NULL)
        {
            ret = 2;
        }
    }
    
    if (ret == 0)
    {
        // Makes building the model faster if it is done rows by row.
        set_add_rowmode(lp, TRUE);
        
        // Construct first row (120 x + 210 y <= 15000).
        j   = 0;
        
        // First column.
        colno[j]    = 1;
        row[j++]    = 120;
        
        // Second column.
        colno[j]    = 2;
        row[j++]    = 210;
        
        // Add the row to lpsolve.
        if (!add_constraintex(lp, j, row, colno, LE, 15000))
        {
            ret = 3;
        }
    }
    
    // Construct second row (110x + 30y <= 4000).
    if (ret == 0)
    {
        j = 0;
        
        // First column.
        colno[j]    = 1;
        row[j++]    = 110;
        
        // Second column.
        colno[j]    = 2;
        row[j++]    = 30;
        
        // Add the row to lpsolve.
        if (!add_constraintex(lp, j, row, colno, LE, 4000))
        {
            ret = 3;
        }
    }
    
    // Construct third row (x + y <= 75).
    if (ret == 0)
    {
        j = 0;
        
        // First column.
        colno[j]    = 1;
        row[j++]    = 1;
        
        // Second column.
        colno[j]    = 2;
        row[j++]    = 1;
        
        // Add the row the lpsolve.
        if (!add_constraintex(lp, j, row, colno, LE, 75))
        {
            ret = 3;
        }
    }
    
    // Set the objective function (143x + 60y).
    if (ret == 0)
    {
        // Row mode should be truned off again when done building the model.
        set_add_rowmode(lp, FALSE);
        
        // Set the objective function (143x + 60y).
        j = 0;
        
        // First column.
        colno[j]    = 1;
        row[j++]    = 143;
        
        // Second column.
        colno[j]    = 2;
        row[j++]    = 60;
        
        // Set the objective in lpsolve.
        if (!set_obj_fnex(lp, j, row, colno))
        {
            ret = 4;
        }
    }
    
    if (ret == 0)
    {
        // Set the object direction to maximize.
        set_maxim(lp);

        // Just out of curioucity, now show the model in lp format on screen.
        // This only works if this is a console application. If not, use
        // wirte_lp and a file name.
        write_LP(lp, stdout);
        
        // I only want to see important messages on screen while solving.
        set_verbose(lp, IMPORTANT);
        
        // Now let lpsolve calculate a solution.
        ret = solve(lp);
        
        if (ret == OPTIMAL)
        {
            ret = 0;
        }
        else
        {
            ret = 5;
        }
    }
    
    // A solution is calculated, now lets get some results.
    if (ret == 0)
    {
        // Objective value.
        cout<<"Objective value: "<<get_objective(lp)<<endl;
        
        // Variable values.
        get_variables(lp, row);
        
        for (j = 0; j < Ncol; j++)
        {
            cout<<get_col_name(lp, j+1)<<"="<<row[j]<<endl;
        }
        
        // We are done now.
    }
    
    // Free allocated memory.
    if (row != NULL)
    {
        free(row);
    }
    
    if (colno != NULL)
    {
        free(colno);
    }
    
    // Clean up such that all used memory by lpsolve is freed.
    if (lp != NULL)
    {
        delete_lp(lp);
    }
    
    return ret;
    
}

void readPtsFile(string filename){
    ifstream infile;
    infile.open(filename);
    if(!infile)
    {
        cout << "fail to open file " << filename << endl;
    }
    string str;
    while(getline(infile,str))   //按行读取,遇到换行符结束
    {
        cout<<str<<endl;
    }
    infile.close();
}

void readPtsFileCVersion(char* filename){
    FILE *fp = std::fopen(filename,"r+");
//    for (int k = 0; k<3; k++) fscanf(f, "%*[^\n]%*c");//跳过前三行，用到了正则表达式
//    Mat1f temp(2 * 68, 1);//点云存储在一个float类型的向量中。
//    for (int j = 0; j<68; j++){
//        fscanf(f, "%f", &temp(2 * j, 0));//依次读取每一行中的float数值，并存储下来
//        fscanf(f, "%f", &temp(2 * j + 1, 0));
//    }
    if(NULL == fp){
        printf("Open Failed.\n");
    }
    char StrLine[1024];             //每行最大读取的字符数
    while (!feof(fp))
    {
        fgets(StrLine,1024,fp);  //读取一行
        printf("%s\n", StrLine); //输出
    }
}

int countLines(string filename){
    ifstream ReadFile;
    int n=0;
    string tmp;
    //ReadFile.open(filename,ios::in);//ios::in 表示以只读的方式读取文件
    ReadFile.open(filename);
    if(ReadFile.fail())//文件打开失败:返回0
    {
        cout << "fail to open in count line" << endl;
        ReadFile.close();
        return 0;
    }
    else//文件存在
    {
        while(getline(ReadFile,tmp))
        {
            n++;
        }
        cout << n << endl;
        ReadFile.close();
        return n;
    }
}

int ** creatPointCloudArray(int rows){
    int **pointcloud = new int*[rows];
    for (int i = 0; i < rows; i++){
        pointcloud[i] = new int[3];
    }
    return pointcloud;
}

void deletePointCloudArray(int **pointcloud, int rows){
    for (int i = 0; i < rows; i++){
        delete []pointcloud[i];
        pointcloud[i] = NULL;
    }
    delete []pointcloud;
    pointcloud = NULL;
}

vector<vector<double>> createPointCloudArray(int rows){
    vector<vector<double>> a(rows,vector<double>(3));
    return a;
}

void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;
    while(std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2-pos1));
        
        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if(pos1 != s.length())
        v.push_back(s.substr(pos1));
}

vector<vector<double>> creatPointCloudArrayFromFile(string filename){
    
    vector<vector<double>> pointclouds;
    
    ifstream infile;
    infile.open(filename);
    if(!infile)
    {
        cout << "fail to open file " << filename << endl;
    }
    
    string str;
    std::string::size_type pos1, pos2;
    const string space = " ";
    
    while(getline(infile,str))   //按行读取,遇到换行符结束
    {
        //cout<<str<<endl;
        
        vector<double> temp;
        
        pos2 = str.find(space);
        pos1 = 0;
        while(std::string::npos != pos2)
        {
            temp.push_back(atof((str.substr(pos1, pos2-pos1)).c_str()));
            
            pos1 = pos2 + space.size();
            pos2 = str.find(space, pos1);
        }
        if(pos1 != str.length())
            temp.push_back(atof((str.substr(pos1)).c_str()));
        
        pointclouds.push_back(temp);
    }
    infile.close();
    return pointclouds;
}

void printVector(vector<vector<double>> pointclouds){
    long count = pointclouds.size();
    for (int i = 0; i < count ; i++)
    {
        cout << pointclouds[i][0] << " " << pointclouds[i][1] << " " << pointclouds[i][2] << endl;
    }
}
