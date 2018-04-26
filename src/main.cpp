/*
统计输出的次数
*/
#include <iostream>
#include <memory.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <functional>
#include <ctime>
#include <map>
#include <set>

using namespace std;

//----------------------------------STBF本身的参数

double hash_A = (sqrt(5) - 1) / 2;
double hash_B = (sqrt(7) - 2) / 2;
double hash_C = sqrt(3) / 2;


#define a 300			//周期数量
#define b 2000		//每个周期包含的元素个数
#define m 2000 	//一个STBF包含的格子个数(for hot)
#define w1 192000	//底层counter数目
#define hash_print_max 0xFFFF	//hash_print最大值
#define hash_print_1_max 0xFF	//另一个指纹最大值
#define int_max 0xFFFFFFFF	//unsigned int的最大值
#define g 110				//假设一个threshold g=bit_of_ID / bit_of_Raptor
#define bit_of_ID 32	//ID的长度 4 * 8
#define str_max_len 4	//读入字符串长度
#define STR_MAX_LEN_INPUT 8

//---------------------------cold filter需要的参数


#define th1 15	//底层threshold
#define hash_number	3//hash函数个数


//------------------------------------
//------------------------------------



//4个byte
class STBF {
public:
	bool raptor;
	unsigned short hash_print;
	unsigned char hash_print_1;
	bool flag;
	STBF() {
		raptor = 0;
		hash_print = 0;
		flag = 0;
		hash_print_1 = 0;
	};
};


class ID {
public:
	char x[str_max_len] = { 0 };
	bool if_equal(ID n) {
		bool flag = 1;
		for (int i = 0; i < str_max_len;i++) {
			if (x[i] != n.x[i]) {
				flag = 0;
				break;
			}
		}
		return flag;
	}
};

class ID_input {
public:
	char x[STR_MAX_LEN_INPUT] = { 0 };
};

bool operator < (ID an, ID bn) {
	for (int i = 0;i < str_max_len;i++) {
		if (bn.x[i] < an.x[i]) {
			return true;
		}
		else if (bn.x[i] > an.x[i]) {
			return false;
		}
	}
	return false;
}


ID_input data_of_id[a * b] = { 0 };
bool hash_print_memory[hash_print_max] = { 0 };	//用于记录某个指纹是否已经尝试解码
unsigned char L1[w1] = { 0 };				//底层counter
STBF s[a * m];		//STBF for hot
multimap<ID, int> id_map;
map<ID, int> all_id_map;
set<ID> id_set;
//字符串hash
uint32_t murmur3_32(const char* key, size_t len, uint32_t seed) {
	uint32_t h = seed;
	if (len > 3) {
		const uint32_t* key_x4 = (const uint32_t*)key;
		size_t i = len >> 2;
		do {
			uint32_t k = *key_x4++;
			k *= 0xcc9e2d51;
			k = (k << 15) | (k >> 17);
			k *= 0x1b873593;
			h ^= k;
			h = (h << 13) | (h >> 19);
			h = (h * 5) + 0xe6546b64;
		} while (--i);
		key = (const char*)key_x4;
	}
	if (len & 3) {
		size_t i = len & 3;
		uint32_t k = 0;
		key = &key[i - 1];
		do {
			k <<= 8;
			k |= *key--;
		} while (--i);
		k *= 0xcc9e2d51;
		k = (k << 15) | (k >> 17);
		k *= 0x1b873593;
		h ^= k;
	}
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

//定义三个hash函数
int hash_1(int hash_range, unsigned int ID) {
	double x = ID * hash_A;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
};

int hash_2(int hash_range, unsigned int ID) {
	double x = ID * hash_B;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
}

int hash_3(int hash_range, unsigned int ID) {
	double x = ID * hash_C;
	double y = x - (unsigned int)x;
	int z = (int)(hash_range * y);
	return z;
}

//指纹
unsigned short hash_o(unsigned int ID) {
	double x = ID * hash_A;
	double y = x - (unsigned int)x;
	unsigned short z = (unsigned short)(hash_print_max * y);
	return z;
}

//另一个指纹
unsigned char hash_o_1(unsigned int ID) {
	double x = ID * hash_B;
	double y = x - (unsigned int)x;
	unsigned char z = (unsigned char)(hash_print_1_max * y);
	return z;
}

//如果是hot元素返回1，cold元素返回0
bool cold_filter(unsigned int ID) {
	unsigned char V1 = th1;
	int h_L1[hash_number];
	h_L1[0] = hash_1(w1, ID);
	h_L1[1] = hash_2(w1, ID);
	h_L1[2] = hash_3(w1, ID);
	for (int j = 0; j < hash_number; j++) {
		if (L1[h_L1[j]] < V1) {
			V1 = L1[h_L1[j]];
		}
	}
	if (V1 < th1) {
		for (int j = 0; j < hash_number; j++) {
			if (L1[h_L1[j]] == V1) {
				L1[h_L1[j]] = L1[h_L1[j]] + 1;
			}
		}
		return 0;
	}
	else {
		return 1;
	}
}

//编码
bool encode(ID myid, int seed) {
	srand(seed);
	bool coeff_a[bit_of_ID];
	for (int i = 0; i < bit_of_ID;i++) {
		coeff_a[i] = rand() % 2;
	}



	bool ID_binary[bit_of_ID];
	for (int i = 0; i < str_max_len; i++) {
		for (int j = 0; j < 8;j++) {
			ID_binary[i * 8 + j] = bool(myid.x[i] & (1 << (8 - 1 - j)));
		}
	}
	bool Raptor = 0;
	for (int i = 0; i < bit_of_ID; i++) {
		Raptor = Raptor ^ (coeff_a[i] & ID_binary[i]);
	}

	return Raptor;
}

//Gauss消元
ID gauss(int _cnt, bool long_coeff[][bit_of_ID], bool *long_raptor) {

	//高斯若儿当消元

	for (int i = 0; i < bit_of_ID; i++) {
		if (long_coeff[i][i] == 0) {
			for (int j = i + 1; j < _cnt; j++) {
				if (long_coeff[j][i] == 1) {
					for (int k = 0; k < bit_of_ID; k++) {
						long_coeff[i][k] = long_coeff[i][k] ^ long_coeff[j][k];
					}
					long_raptor[i] = long_raptor[i] ^ long_raptor[j];
					break;
				}
			}
		}

		for (int j = i + 1; j < _cnt; j++) {
			if (long_coeff[j][i] == 1) {
				for (int k = 0; k < bit_of_ID; k++) {
					long_coeff[j][k] = long_coeff[j][k] ^ long_coeff[i][k];
				}
				long_raptor[j] = long_raptor[j] ^ long_raptor[i];
			}

		}

	}



	for (int i = bit_of_ID - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			if (long_coeff[j][i] == 1) {
				long_coeff[j][i] = long_coeff[j][i] ^ long_coeff[i][i];
				long_raptor[j] = long_raptor[j] ^ long_raptor[i];
			}
		}
	}



	ID myid;


	for (int i = 0; i < str_max_len; i++) {
		for (int j = 0; j < 8;j++) {
			if (long_raptor[i * 8 + j]) {
				myid.x[i] = myid.x[i] | (1 << (8 - 1 - j));
			}
		}
	}

	return myid;
}



//构造方程的系数
void parameter(int &_cnt, int *cnt, bool *raptor, bool  long_coeff[][bit_of_ID], bool * long_raptor) {
	for (int i = 0; i < _cnt; i++) {
		srand(cnt[i]);
		for (int j = 0; j < bit_of_ID;j++) {
			long_coeff[i][j] = rand() % 2;
		}
	}

	for (int i = 0; i < _cnt; i++) {
		long_raptor[i] = raptor[cnt[i]];
	}

}

//记录数据到STBF
void record(int _m, int i, STBF *s, unsigned int id_num, ID myid) {

	int h1 = hash_1(_m, id_num);
	int h2 = hash_2(_m, id_num);
	int h3 = hash_3(_m, id_num);
	unsigned short ho = hash_o(id_num);
	unsigned char ho_1 = hash_o_1(id_num);
	bool myraptor = encode(myid, i);
	//没有记录
	if (s[i * _m + h1].flag == 0 && s[i * _m + h1].raptor == 0 && s[i * _m + h1].hash_print == 0 && s[i * _m + h1].hash_print_1 == 0) {
		s[i * _m + h1].flag = 1;
		s[i * _m + h1].raptor = myraptor;
		s[i * _m + h1].hash_print = ho;
		s[i * _m + h1].hash_print_1 = ho_1;
	}
	//发生碰撞
	else if (s[i * _m + h1].flag == 1 && (s[i * _m + h1].raptor != myraptor || s[i * _m + h1].hash_print != ho || s[i * _m + h1].hash_print_1 != ho_1)) {
		s[i * _m + h1].flag = 0;
		s[i * _m + h1].raptor = 1;
		s[i * _m + h1].hash_print = hash_print_max;
		s[i * _m + h1].hash_print_1 = hash_print_1_max;
	}

	if (s[i * _m + h2].flag == 0 && s[i * _m + h2].raptor == 0 && s[i * _m + h2].hash_print == 0 && s[i * _m + h2].hash_print_1 == 0) {
		s[i * _m + h2].flag = 1;
		s[i * _m + h2].raptor = myraptor;
		s[i * _m + h2].hash_print = ho;
		s[i * _m + h2].hash_print_1 = ho_1;
	}
	else if (s[i * _m + h2].flag == 1 && (s[i * _m + h2].raptor != myraptor || s[i * _m + h2].hash_print != ho || s[i * _m + h2].hash_print_1 != ho_1)) {
		s[i * _m + h2].flag = 0;
		s[i * _m + h2].raptor = 1;
		s[i * _m + h2].hash_print = hash_print_max;
		s[i * _m + h2].hash_print_1 = hash_print_1_max;
	}

	if (s[i * _m + h3].flag == 0 && s[i * _m + h3].raptor == 0 && s[i * _m + h3].hash_print == 0 && s[i * _m + h3].hash_print_1 == 0) {
		s[i * _m + h3].flag = 1;
		s[i * _m + h3].raptor = myraptor;
		s[i * _m + h3].hash_print = ho;
		s[i * _m + h3].hash_print_1 = ho_1;
	}
	else if (s[i * _m + h3].flag == 1 && (s[i * _m + h3].raptor != myraptor || s[i * _m + h3].hash_print != ho || s[i * _m + h3].hash_print_1 != ho_1)) {
		s[i * _m + h3].flag = 0;
		s[i * _m + h3].raptor = 1;
		s[i * _m + h3].hash_print = hash_print_max;
		s[i * _m + h3].hash_print_1 = hash_print_1_max;
	}
}


//对每一列中的元素分组
void form_group_col(vector<vector<int>> &G, int _m, int col_i, STBF *s) {
	for (int i = 0;i < a;i++) {
		if (s[i * _m + col_i].flag == 1 && hash_print_memory[s[i * _m + col_i].hash_print] == 0) {
			if (G.empty()) {
				vector<int> tmp;
				int tmp_p;
				tmp_p = i;
				tmp.push_back(tmp_p);
				G.push_back(tmp);
			}
			else {
				bool flag = 0;
				vector<vector<int>>::iterator iter;
				for (iter = G.begin(); iter != G.end(); iter++) {
					if (s[i * _m + col_i].hash_print == s[*iter->begin() * _m + col_i].hash_print && s[i * _m + col_i].hash_print_1 == s[*iter->begin() * _m + col_i].hash_print_1) {
						iter->push_back(i);
						flag = 1;
						break;
					}
				}
				if (flag == 0) {
					vector<int> tmp;
					int tmp_p;
					tmp_p = i;
					tmp.push_back(tmp_p);
					G.push_back(tmp);
				}
			}
		}
	}
}

//寻找该列有几个碰撞元素
int find_colli(int _m, int col_i, STBF *s) {
	int cnt = 0;
	for (int i = 0;i < a;i++) {
		if (s[i * _m + col_i].flag == 0 && s[i * _m + col_i].hash_print == hash_print_max && s[i * _m + col_i].raptor == 1) {
			cnt++;
		}
	}
	return cnt;
}

//记录一共有几行出现_cnt，以及cnt[]存储行标。cnt[i]存储顺序地出现地行数，raptor[i]存储实际i行里地raptor。
void record_rows(int &_cnt, int *cnt, int col_i, bool *raptor, vector<vector<int>>::iterator iter, STBF *s) {
	bool memory_row[a] = { 0 }; //记录某一行是否有该元素
	vector<int>::iterator iter1;
	for (iter1 = iter->begin(); iter1 != iter->end(); iter1++) {
		if (memory_row[*iter1] == false) {
			memory_row[*iter1] = true;
			raptor[*iter1] = s[*iter1 * m + col_i].raptor;
		}
	}
	for (int i = 0; i < a; i++) {
		if (memory_row[i] == 1) {
			cnt[_cnt] = i;
			_cnt++;
		}
	}
}



int main() {
	int out_cnt = 0;
	//读入数据并存储
	fstream fin("C:\\Users/Corta/Desktop/1806/persistent_item_data/data_subset/new_zipf/030.dat", ios::in | ios::binary);
	//ofstream draw_out;
	//draw_out.open("C:\\Users/Corta/Desktop/1806/persistent_item_data/plot/23.dat", ios::app);
	//ofstream draw_out1;
	//draw_out1.open("C:\\Users/Corta/Desktop/1806/persistent_item_data/plot/15.dat", ios::app);
	unsigned int i = 0;
	for (i = 0; i < a * b; i++) {
		fin.read((char *)(&data_of_id[i]), sizeof(data_of_id[i]));
	}
	fin.close();
	int in_cnt = 0;
	unsigned int j = 0;
	for (i = 0;i < a;i++) {
		id_set.clear();
		for (j = 0;j < b;j++) {
			ID input_my_id;
			memcpy(input_my_id.x, data_of_id[i * b + j].x, sizeof(input_my_id.x));
			unsigned int _data = murmur3_32(input_my_id.x, str_max_len, 1);
			if (cold_filter(_data)) {
				record(m, i, s, _data, input_my_id);
				in_cnt++;
			}
			if (id_set.find(input_my_id) == id_set.end()) {
				id_set.insert(input_my_id);
				all_id_map[input_my_id] += 1;
			}
		}
	}
	cout << in_cnt << endl;

	int BF_cnt_0 = 0;
	int BF_cnt_1 = 0;
	int BF_cnt_colli = 0;
	for (i = 0; i < a; i++) {
		for (j = 0; j < m;j++) {
			if (s[i * m + j].flag == 0 && s[i * m + j].hash_print != hash_print_max) {
				BF_cnt_0++;
			}
			else {
				BF_cnt_1++;
			}
			if (s[i * m + j].flag == 0 && s[i * m + j].hash_print == hash_print_max && s[i * m + j].hash_print_1 == hash_print_1_max) {
				BF_cnt_colli++;
			}
		}
	}

	cout << "空  " << BF_cnt_0 << endl;
	cout << "有效" << BF_cnt_1 << endl;
	cout << "撞  " << BF_cnt_colli << endl;
	system("pause");
	cout << all_id_map.size() << endl;
	system("pause");
	map<ID, int>::iterator it;
	int acc_cnt = 0;
	for (it = all_id_map.begin();it != all_id_map.end();it++) {
		if (it->second > g) {
			acc_cnt++;
		}
	}
	int col_i = 0;
	for (col_i = 0; col_i < m; col_i++) {
		int num_coli = find_colli(m, col_i, s);
		//对hot每一列分组，寻找是否可解
		vector<vector<int>> G;	//G记录相同指纹的行标
		form_group_col(G, m, col_i, s);
		vector<vector<int>>::iterator iter;
		for (iter = G.begin();iter != G.end();iter++) {
			if (iter->size()  > g) {
				//hash_print_memory[s[*iter->begin() * m + col_i].hash_print][s[*iter->begin() * m + col_i].hash_print_1] = 1;	//可解时，记录已被解
				int _cnt = 0;
				int cnt[a] = { 0 };
				bool raptor[a] = { 0 };
				record_rows(_cnt, cnt, col_i, raptor, iter, s);
				int coli[2] = { 0 };
				coli[0] = -1;
				coli[1] = -1;

				//寻找碰撞元素
				for (i = 0; i < a; i++) {
					if (s[i * m + col_i].flag == 0 && s[i * m + col_i].hash_print == hash_print_max && s[i * m + col_i].hash_print_1 == hash_print_1_max) {
						if (coli[0] == -1 || coli[1] == -1) {
							for (j = 0; j < m;j++) {
								if (s[*iter->begin() * m + col_i].hash_print == s[i * m + j].hash_print && s[*iter->begin() * m + col_i].hash_print_1 == s[i * m + j].hash_print_1) {
									cnt[_cnt] = i;
									raptor[i] = s[i * m + j].raptor;
									_cnt++;

									if (coli[0] == -1) {
										coli[0] = j;
									}
									else {
										if (coli[0] != j && coli[1] == -1) {
											coli[1] = j;
										}
									}
									break;
								}
							}
						}
						else {
							for (int coli_i = 0; coli_i < 2;coli_i++) {
								if (s[*iter->begin() * m + col_i].hash_print == s[i * m + coli[coli_i]].hash_print && s[*iter->begin() * m + col_i].hash_print_1 == s[i * m + coli[coli_i]].hash_print_1) {
									cnt[_cnt] = i;
									raptor[i] = s[i * m + coli[coli_i]].raptor;
									_cnt++;
									break;
								}
							}
						}
					}
				}

				//hot中该个数足够解出线性方程
				if (_cnt > g) {
					//存储系数
					bool *long_raptor = new bool[_cnt]();


					bool(*long_coeff)[bit_of_ID] = new bool[_cnt][bit_of_ID]();
					parameter(_cnt, cnt, raptor, long_coeff, long_raptor);

					//高斯若儿当消元

					ID myid = gauss(_cnt, long_coeff, long_raptor);
					unsigned int id_num = murmur3_32(myid.x, str_max_len, 1);
					//检查ID是否正确
					bool test_flag = 1;
					unsigned short test_ho = hash_o(id_num);
					unsigned char test_ho_1 = hash_o_1(id_num);
					int h1 = hash_1(m, id_num);
					int h2 = hash_2(m, id_num);
					int h3 = hash_3(m, id_num);
					if (test_ho == s[*iter->begin() * m + col_i].hash_print && test_ho_1 == s[*iter->begin() * m + col_i].hash_print_1) {
						if (col_i != h1 && col_i != h2 && col_i != h3) {
							test_flag = 0;
						}
					}
					else {
						test_flag = 0;
					}
					if (test_flag) {
						out_cnt++;
						id_map.insert(make_pair(myid, _cnt));
						hash_print_memory[s[*iter->begin() * m + col_i].hash_print] = 1;
					}
					delete[] long_raptor;
					delete[] long_coeff;
				}
			}
		}
	}

	for (multimap<ID, int>::iterator iter = id_map.begin(); iter != id_map.end(); iter++) {

		double x1 = all_id_map[iter->first];
		double x2 = iter->second;
		//draw_out << "PIE+CF," << (double) (4000 - m) / (double) 4000 << "," << x1 - x2 << endl;
		//draw_out1 << "PIE+CF," << (double)(4000 - m) / (double)4000 << "," << 1 - (x2/x1) << endl;

	}
	double recall_rate = (double)out_cnt / (double)acc_cnt;
	cout << out_cnt << endl;
	system("pause");
	//draw_out << "$\mathcal{T}_1 =$" << th1 << "," << (double)m * a * 4 / 1000000 << "," << recall_rate << endl;
	//draw_out.close();
	//draw_out << "PIE+CF," << (double)(4000 - m) / (double)4000 << "," << recall_rate << endl;
	//draw_out.close();
	//draw_out1.close();

	return 0;
}
