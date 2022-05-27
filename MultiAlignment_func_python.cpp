#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <math.h>
using namespace std;

//g++ -fPIC MultiAlignment_func_python.cpp -o MultiAlignment_func_python.so -shared -I/home/limengtian/miniconda3/envs/python3/include/python3.8

struct Result_merge{
    int score;
    int start_i;
    int end_i;
    int start_j;
    int end_j;
    char *seq;
    int length;
};

struct Result_merge* MultiOptimal_merge_alignment (char *seq_a, char *seq_b, 
                                    double plus=1, double mu=2, double delta=10, int multi = 3) {
    int i_max=0;
    int j_max=0;

    int L_a=strlen(seq_a);
    int L_b=strlen(seq_b);
    int *H[L_a+1];
    int *Label[L_a+1];
    
    for(int i=0;i<=L_a;i++){
        H[i] = new int[L_b+1]; 
        Label[i] = new int[L_b+1]; 
        for(int j=0;j<=L_b;j++){
            H[i][j]=0;
			Label[i][j]=0;
        }
    }

    double temp[4];
    int *trace[L_a+1]; 

    for(int i=0;i<=L_a;i++){
        trace[i] = new int[L_b+1]; 
        for(int j=0;j<=L_b;j++){
            trace[i][j]=3;
        }
    }

    for(int i = 1; i <= L_a; i++){
        for(int j = 1; j <= L_b; j++){
            if(seq_a[i-1] == seq_b[j-1])  temp[0] = H[i-1][j-1] + plus;
            else  temp[0] = H[i-1][j-1] - mu;
            temp[1] = H[i-1][j] - delta;
            temp[2] = H[i][j-1] - delta;
            temp[3] = 0;

            double max = temp[0];
            trace[i][j] = 0;
            
            for(int k = 1; k < 4; ++k){
                if(temp[k] > max){
                    max = temp[k];
                    trace[i][j] = k;
                }
            }
            H[i][j] = max;
        }
    }

    double H_max = 0;

    for(int i = 1; i <= L_a; i++){
        for(int j = 1; j <= L_b; j++){
                if(H[i][j] > H_max){
                    H_max = H[i][j];
                    i_max = i;
                    j_max = j;
                }
        }
    }
	
    int start_i,start_j;
    int score = H[i_max][j_max];
    int current_i = i_max, current_j = j_max;
    int end_i = i_max, end_j = j_max;
    string traceback_a="";
    string traceback_b="";

    int ***MH = new int** [multi];// Statement: MH[multi][L_a+1][L_b+1];
	
	for(int k = 0; k < multi; k++){
		MH[k] = new int* [L_a+1];
		for(int i=0;i<=L_a;i++){
			MH[k][i] = new int [L_b+1];
		}
	}
	Result_merge* result_multi = new Result_merge[multi];
	
	for(int k = 0; k < multi; k++){
		//Initialize MH[k]
		for(int i = 0; i <= L_a; i++){
			for(int j = 0; j <= L_b; j++){
				if (k==0){
					MH[k][i][j]=H[i][j];	
				}
				else{
					MH[k][i][j]=MH[k-1][i][j];	
				}
			}
		}
		
		traceback_a="";
		traceback_b="";
		while(trace[current_i][current_j]!=3){
			MH[k][current_i][current_j] = 0;
			Label[current_i][current_j] = -1;
			switch(trace[current_i][current_j]){
				case 0:
					traceback_a += seq_a[current_i-1];
					traceback_b += seq_b[current_j-1];
					current_i = current_i-1;
					current_j = current_j-1;
					break;
				case 1:
					traceback_b += "N";
					traceback_a += seq_a[current_i-1];
					current_i = current_i-1;
					break;
				case 2:
					traceback_a += "N";
					traceback_b += seq_b[current_j-1];
					current_j = current_j-1;
					break;
			}
			//cout << "Alignment" << current_i <<" "<< current_j<< endl;
		}

		start_i = current_i;
		start_j = current_j;
		//cout<<"start_i:"<<start_i<<" "<<"start_j:"<<" "<<start_j<<endl;

		int length = traceback_a.size();
		//cout<<"length:"<<length<<endl;
		string s = "";
		for (int i = length-1; i >= 0; --i){
			if (traceback_a[i]=='N')  s+=traceback_b[i];
			else  s+=traceback_a[i];
		}

		length = s.size();
		char *seq = new char[length];
		for (int i = 0; i < length; ++i){
			seq[i]=s[i];
		}
		/*/
		result_multi[k] = {score, start_i, end_i, start_j, end_j, seq, length};
		/*/
		result_multi[k].score=score;
		result_multi[k].start_i=start_i+1;
		result_multi[k].end_i=end_i;
		result_multi[k].start_j=start_j+1;
		result_multi[k].end_j=end_j;
		result_multi[k].seq=seq;
		result_multi[k].length=length;
		
		//if (k==multi-1)	break;
		if (start_i>0)		temp[1] = MH[k][start_i-1][start_j]-delta;
		else temp[1] = 0;
		if (start_j>0)	temp[2] = MH[k][start_i][start_j-1]-delta;
		else temp[2] = 0;
		temp[3] = 0;
		double tmax = temp[1];
		trace[start_i][start_j] = 3;
		for(int k = 1; k < 4; ++k){
			if(temp[k] > tmax){
				tmax = temp[k];
				trace[start_i][start_j] = k;
			}
		}
		MH[k][start_i][start_j] = tmax;
		
		//cout << "restart" << start_i <<" "<< start_j<<" "<<tmax<< endl;
		for(int i = start_i+1; i <= L_a; i++){
			for(int j = start_j+1; j <= L_b; j++){
				if(Label[i][j]==0 and seq_a[i-1] == seq_b[j-1])  temp[0] =MH[k][i-1][j-1] + plus;
				else  temp[0] =MH[k][i-1][j-1] - mu;
				temp[1] =MH[k][i-1][j] - delta;
				temp[2] =MH[k][i][j-1] - delta;
				temp[3] = 0;

				tmax = temp[0];
				trace[i][j] = 0;
				
				for(int k = 1; k < 4; ++k){
					if(temp[k] > tmax){
						tmax = temp[k];
						trace[i][j] = k;
					}
				}
				MH[k][i][j]= tmax;
				//if (k==2)	cout << i << " " << j << " " << H[i][j] << " " <<MH[k][i][j] << endl;
				
				/*/if (MH[k][i][j]==H[i][j] and j-start_j>i-start_i)  {
					//cout<< "i"<< i<<" "<<"j" <<j<<endl;
					break;
				}/*/
			}
		}
		
		H_max = 0;
		i_max = 1;
		j_max = 1;
		for(int i = 0; i <= L_a; i++){
			for(int j = 0; j <= L_b; j++){
					if(MH[k][i][j] >H_max){
						H_max = MH[k][i][j];
						i_max = i;
						j_max = j;
					}
			}
		}
		score = MH[k][i_max][j_max];
		current_i = i_max, current_j = j_max;
		end_i = i_max, end_j = j_max;
		//cout<< end_i <<" "<< end_j << " " << score << endl;
	}
	
	for (int k = 0; k < multi; k++){
		for (int i = 0; i < L_a+1; i++){
			delete[]MH[k][i];
		}
		delete[]MH[k];
	}
	delete[]MH;
	for(int i = 0; i <= L_a; ++i){
        delete [] H[i]; delete [] trace[i];delete[] Label[i];
    }

    return result_multi;
}

static PyObject * Conver_MultiOptimal_merge_alignment(PyObject *self, PyObject *args)
{
    char *seq_a;
    char *seq_b;
	int multi;
    double plus, mu, delta;

    if(!(PyArg_ParseTuple(args, "ssdddi", &seq_a, &seq_b,&plus,&mu,&delta,&multi))){
        return NULL;
    }

    Result_merge* result_list = new Result_merge[multi];
    result_list = MultiOptimal_merge_alignment(seq_a, seq_b,plus,mu,delta,multi);
    
    PyObject *Lst = PyList_New(multi);

    for (int i = 0; i < multi; ++i){

        Result_merge result;
        result = result_list[i];
        int score=result.score;
        int start_i=result.start_i;
        int end_i=result.end_i;
        int start_j=result.start_j;
        int end_j=result.end_j;
        char *seq=result.seq;
        int length=result.length;

        PyObject *lst = Py_BuildValue("[i,s#,i,i,i,i]", score, seq, length, start_i, end_i, start_j, end_j);
        PyList_SET_ITEM(Lst, i, lst);
    }

    return Lst;
}

// 定义模组方法
static PyMethodDef module_methods[] = {
    //{"banded_overlap_alignment", Conver_banded_overlap_alignment, METH_VARARGS},
    {"MultiOptimal_merge_alignment", Conver_MultiOptimal_merge_alignment, METH_VARARGS},
    {NULL, NULL},
};

// 定义模组
static struct PyModuleDef MultiAlignment_func_python = {
    PyModuleDef_HEAD_INIT,"MultiAlignment_func_python","",-1,module_methods

};

extern "C"

// 模组初始化
PyMODINIT_FUNC PyInit_MultiAlignment_func_python()
{
     return PyModule_Create(&MultiAlignment_func_python);
}