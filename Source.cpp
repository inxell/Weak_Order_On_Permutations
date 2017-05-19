#include<string>
#include<vector>
#include<algorithm>
#include<iostream>
#include<map>
#include<set>
#include<queue>
#include<stack>
#include<istream>
#include<sstream>
#include<random>
#include<unordered_set>
#include<unordered_map>
#include <iomanip>

using namespace std;

void print_vec(const vector<int>& vec) {
	for (int i : vec) {
		cout << setw(2) << i << ' ';
	}
	cout << endl;
}

void print_vec(const vector<short>& vec) {
	for (short i : vec) {
		cout << setw(2) << i << ' ';
	}
	cout << endl;
}

void print_vec(const vector<long double>& vec) {
	for (long double i : vec) {
		cout << setw(5) << i << ' ';
	}
}

vector<int> invs_to_permut(vector<short> & invs) {
	/*
	Convert vector of inversions to permutation.
	For example, invs = [0, 1, 0, 3, 2] maps to 
	the permutation permut = [4, 2, 5, 1, 3].
	Note that invs[i] denotes #{ k < i : permut[k] > permut[i]}.
	*/
	vector<int> res{ 1 };
	for (short i = 1; i < invs.size(); ++i) {
		short k = i + 1 - invs[i];
		for (int& j : res) {
			if (j >= k) ++j;
		}
		res.push_back(k);
	}
	return res;
}

class Permutation {
public:
	vector<short> invs;
	/*
	invs is a vector with contains the information of the pairs of inversions 
	which determines the true value of the permutation. Specifically, if
	(i, j) is a pair of inversions for the permutation permut, i.e., i < j and
	permut[i] > permut[j], we have invs[i-1] & (1 << j-1) == 1.
	*/
	short len;
	// len is the size of the permutation.
	unsigned cnt_invs = 0;
	// cnt_invs is the number of inversions of the permutation.

	Permutation(short n) : len(n) {
		// construct the identity permutation of size n.
		if (n > 16) {
			cout << "Size Overflow!" << endl;
		}
		invs = vector<short>(n, 0);
		cnt_invs = 0;
	}

	Permutation(const vector<short>& vec) : invs(vec) {
		// construct the permutation with the given vector denotes the pairs of inversions.
		len = vec.size();
		unsigned tmp = 0;
		for (int i = 0; i < len; ++i) {
			short k = vec[i];
			for (int j = i + 1; j < len; ++j) {
				if (k & (1 << j)) {
					++tmp;
				}
			}
		}
		cnt_invs = tmp;
	}

	Permutation(const vector<int>& vec) {
		// construct the permutation from vector denoting a permutation.
		if (vec.size() > 16) {
			cout << "Size Overflow!" << endl;
			return;
		}
		len = vec.size();
		invs = vector<short>(len, 0);
		unsigned tmp = 0;
		for (int i = 0; i < len; ++i) {
			int cnt = 0;
			for (int j = i + 1; j < len; ++j) {
				if (vec[i] > vec[j]) {
					++cnt;
					invs[i] |= (1 << j);
				}
			}
			if (cnt + i >= len) {
				cout << "Not a Valid Permutation!" << endl;
				return;
			}
			tmp += cnt;
		}
		cnt_invs = tmp;
	}

	vector<int> get_value() const {
		// return the permutation as a vector.
		vector<short> tmp(len, 0);
		for (short i = 0; i < len; ++i) {
			short val = invs[i];
			for (short k = i + 1; k < len; ++k) {
				if (val & (1 << k)) {
					++tmp[k];
				}
			}
		}
		return invs_to_permut(tmp);
	}

	string get_string() const {
		// return the permutation as a string.
		vector<int> tmp = get_value();
		string res;
		for (int i : tmp) {
			res.push_back('0' + i);
		}
		return res;
	}

	Permutation inverse() const {
		// return the inverse of the permutation.
		vector<int> tmp = get_value();
		vector<int> res(len, 0);
		for (int i = 0; i < len; ++i) {
			res[tmp[i] - 1] = i + 1;
		}
		return Permutation(tmp);
	}

	Permutation reversal() const {
		// return the reversal of the permutation.
		vector<int> tmp = get_value();
		reverse(tmp.begin(), tmp.end());
		return Permutation(tmp);
	}


};



bool weak_order_less(const Permutation& a, const Permutation& b) {
	for (short i = 0; i < a.len; ++i) {
		if ((a.invs[i] | b.invs[i]) != b.invs[i]) return false;
	}
	return true;
}

Permutation join(const Permutation& a, const Permutation& b) {
	if (a.len != b.len) {
		cout << "Permutations are of different size!" << endl;
		return Permutation(0);
	}
	short n = a.len;
	vector<short> invs(n, 0);
	for (short i = n - 2; i >= 0; --i) {
		short k = a.invs[i] | b.invs[i];
		for (short j = n - 2; j > i; --j) {
			if (k & (1 << j)) {
				k |= invs[j];
			}
		}
		invs[i] = k;
	}
	return Permutation(invs);
}

Permutation meet(const Permutation& a, const Permutation& b) {
	vector<int> aa = a.get_value(), bb = b.get_value();
	reverse(aa.begin(), aa.end());
	reverse(bb.begin(), bb.end());
	return join(Permutation(aa), Permutation(bb)).reversal();
}


/*
From this below, we code the permutation in a conventional way as a vector of integers.
This way support permutation with larger size. As a trade-off, it is more memory consuming and 
less efficient when checking whether two permutations are comparable in terms of weak order.
*/

bool weak_order_less(const vector<int>& a, const vector<int>& b) {
	for (int i = 0; i < a.size(); ++i) {
		for (int j = i + 1; j < a.size(); ++j) {
			if (a[i] > a[j] && b[i] < b[j]) return false;
		}
	}
	return true;
}

vector<unordered_set<int>> permut_to_pairs(vector<int> p) {
	vector<unordered_set<int>> res(p.size());
	for (int i = 0; i < p.size(); ++i) {
		for (int j = i + 1; j < p.size(); ++j) {
			if (p[i] > p[j]) {
				res[i].insert(j);
			}
		}
	}
	return res;
}

int pairs_to_count_of_invs(const vector<unordered_set<int>>& pairs) {
	int res = 0;
	for (const auto& s : pairs) {
		res += s.size();
	}
	return res;
}

vector<int> reversal(vector<int> p) {
	reverse(p.begin(), p.end());
	return p;
}


vector<int> invs_to_permut(const vector<int>& vec, vector<unordered_set<int>>& pairs) {
	vector<int> res;
	int n = vec.size();
	for (int i = 0; i < vec.size(); ++i) {
		int k = i + 1 - vec[i];
		res.push_back(k);
		for (int j = 0; j < i; ++j) {
			if (res[j] >= k) {
				++res[j];
				pairs[j].insert(i);
			}
		}
	}
	return res;
}


void transitive_closure(vector<unordered_set<int>>& pairs) {
	for (int i = int(pairs.size()) - 2; i >= 0; --i) {
		vector<int> tmp(pairs[i].begin(), pairs[i].end());
		for (int t : tmp) {
			for (int j : pairs[t]) {
				pairs[i].insert(j);
			}
		}
	}
}

vector<int> pairs_to_invs(const vector<unordered_set<int>>& pairs) {
	vector<int> res(pairs.size(), 0);
	for (int i = 0; i < pairs.size(); ++i) {
		for (int k : pairs[i]) {
			++res[k];
		}
	}
	return res;
}

vector<int> pairs_to_permut(const vector<unordered_set<int>>& pairs) {
	vector<int> res, invs = pairs_to_invs(pairs);
	for (int i = 0; i < invs.size(); ++i) {
		if (invs[i] > i) {
			cout << "Invalid Input!" << endl;
			return{};
		}
	}
	for (int i = 0; i < invs.size(); ++i) {
		int k = i + 1 - invs[i];
		res.push_back(k);
		for (int j = 0; j < i; ++j) {
			if (res[j] >= k) {
				++res[j];
			}
		}
	}
	return res;
}



class PermutGen {
public:
	vector<int> gen(int n, vector<unordered_set<int>>& pairs) {
		vector<int> invs;
		pairs = vector<unordered_set<int>>(n);
		for (int j = 1; j <= n; ++j) {
			invs.push_back(dist(mt) * j);
		}
		return invs_to_permut(invs, pairs);
	}

	vector<int> gen(short n) {
		vector<short> invs;
		for (short i = 1; i <= n; ++i) {
			invs.push_back(dist(mt)*i);
		}
		return invs_to_permut(invs);
	}


private:
	random_device rd;
	mt19937 mt = mt19937(rd());
	uniform_real_distribution<> dist = uniform_real_distribution<>(0, 1);
};


// The alterative implementation of permutation and weak order ends here.



vector<int> inverse(const vector<int>& a) {
	// return the inverse of the permutation a.
	vector<int> res(a.size(), 0);
	for (int i = 0; i < a.size(); ++i) {
		res[a[i] - 1] = i + 1;
	}
	return res;
}


struct TreeNode {
	// children are greater than their parents in terms of weak order.
	vector<TreeNode*> children;
	Permutation* val = nullptr;
	TreeNode(Permutation* ptr): val(ptr) {}
};


class Diagram {
public:
	// len denotes the size of the permutation.
	int len = 0;
	// root will be the identity in S_n.
	TreeNode* root = nullptr;
	// graph will point to the private data member visited.
	const map<string, TreeNode*>* graph = nullptr;


	Diagram(short n) {
		/*
		Construct the diagram of S_n under weak order.
		n should be no greater than 16. Practically, n should be less then 8.
		When n == 10, the Diagram will consume 2GB rams and it takes one hour to 
		generate the diagram on my laptop with CPU i5-6300U.
		*/
		if (n <= 0 || n > 16) {
			cout << "Invalid value of n" << endl;
			return;
		}
		len = n;
		auto tmp = new Permutation(n);
		root = new TreeNode(tmp);
		queue<TreeNode*> que;
		visited.insert({ tmp->get_string(), root });
		que.push(root);
		vector<int> pinverse, permut;
		vector<short> invs_vec;
		Permutation* child = nullptr;
		TreeNode* child_node = nullptr;
		while (!que.empty()) {
			TreeNode* curr = que.front();
			que.pop();
			permut = curr-> val -> get_value();
			invs_vec = curr-> val -> invs;
			pinverse = inverse(permut);
			for (int i = 1; i < pinverse.size(); ++i) {
				if (pinverse[i - 1] > pinverse[i]) continue;
				int j = pinverse[i - 1]-1, k = pinverse[i]-1;
				vector<short> new_invs(invs_vec.begin(), invs_vec.end());
				new_invs[j] |= (1 << k);
				child = new Permutation(new_invs);
				string cs = child->get_string();
				if (visited.find(cs) != visited.end()) {
					curr->children.push_back(visited[cs]);
					delete child;
				}
				else {
					child_node = new TreeNode(child);
					curr->children.push_back(child_node);
					visited.insert({ cs, child_node });
					que.push(child_node);
				}

			}
		}

		graph = &visited;
	}

	vector<short> increasing_set_inv_count(vector<string>& permutations) const {
		/*
		return the count of permutations with different number of inversions
		within the increasing set generated by those permutations within the input.
		For example, if the return value is [0, 1, 2, 1], it means that the increasing
		set contains one permutation with 1 inversion, two with 2 inversions and one
		with 3 inversions.
		*/
		unordered_set<TreeNode*> record;
		for (string p : permutations) {
			if (visited.find(p) == visited.end()) {
				cout << "Invalid Permutations" << endl;
				return{};
			}
		}
		for (string p : permutations) {
			TreeNode* curr = visited.at(p); //const member function must use 'at' method
			if (record.find(curr) != record.end()) {
				continue;
			}
			record.insert(curr);
			stack<TreeNode*> stk;
			stk.push(curr);
			while (!stk.empty()) {
				curr = stk.top();
				stk.pop();
				for (TreeNode* child : curr->children) {
					if (record.find(child) == record.end()) {
						stk.push(child);
						record.insert(child);
					}
				}
			}
		}
		int k = len*(len - 1) / 2 + 1;
		vector<short> res(k, 0);
		for (TreeNode* t : record) {
			++res[t->val->cnt_invs];
		}
		return res;
	}

	vector<vector<string>> independent_sets() {
		// generate all the independent sets under weak order.
		// It is impractical on a laptop when n >= 5.
		vector<string> sn;
		for (auto it = visited.begin(); it != visited.end(); ++it) {
			sn.push_back(it->first);
		}
		vector<vector<string>> res;
		for (string p : sn) {
			TreeNode* node_p = visited[p];
			vector<vector<string>> tmp;
			tmp.push_back({ p });
			for (const vector<string>& v : res) {
				bool boo = true;
				for (const string& s : v) {
					TreeNode* ss = visited[s];
					if (weak_order_less(*(node_p->val), *(ss->val)) || weak_order_less(*(ss->val), *(node_p->val))) {
						boo = false;
						break;
					}
				}
				if (boo) {
					vector<string> to_add(v.cbegin(), v.cend());
					to_add.push_back(p);
					tmp.push_back(to_add);
				}
			}
			res.insert(res.end(), tmp.begin(), tmp.end());
		}
		return res;
	}

private:
	map<string, TreeNode*> visited;
};





long double partition_function(int n, long double q) {
	// return the normalizer of the Mallow measure on S_n with parameter q.
	vector<long double> record{ 1 };
	long double power = 1, res = 1;
	for (int i = 1; i < n; ++i) {
		power *= q;
		record.push_back(record.back() + power);
		res *= record.back();
	}
	return res;
}


long double mallows_measure_of_increasing_set(int n, vector<short>& invs, long double q) {
	// Here the invs should be the value returned by the method increasing_set_inv_count of Diagram.
	long double power = 1;
	long double res = 0;
	for (int i = 0; i < invs.size(); ++i) {
		res += invs[i] * power;
		power *= q;
	}
	return res/partition_function(n, q);
}


int main() {
	/*
	// some test code to the class Permutation based on the two different implementation of permutation
	PermutGen per;
	short n = 15;
	int i1, i2, i3, i4;
	for (int ii = 0; ii < 2000; ++ii) {
		vector<int> invs, p1, p2, p3, p4, rp1, rp2;
		vector<unordered_set<int>> pairs1, pairs2, pairs3(n, {}), pairs4(n);
		p1 = per.gen(n, pairs1);
		p2 = per.gen(n, pairs2);
		i1 = pairs_to_count_of_invs(pairs1);
		i2 = pairs_to_count_of_invs(pairs2);
		for (int i = 0; i < n; ++i) {
			for (int j : pairs1[i]) {
				pairs3[i].insert(j);
			}
			for (int j : pairs2[i]) {
				pairs3[i].insert(j);
			}
		}
		transitive_closure(pairs3);
		p3 = pairs_to_permut(pairs3);
		i3 = pairs_to_count_of_invs(pairs3);
		//find the meet of p1 and p2
		rp1 = reversal(p1);
		rp2 = reversal(p2);
		pairs4 = permut_to_pairs(rp1);
		auto tmp = permut_to_pairs(rp2);
		for (int i = 0; i < n; ++i) {
			for (int j : tmp[i]) {
				pairs4[i].insert(j);
			}
		}
		transitive_closure(pairs4);
		i4 = n*(n - 1) / 2 - pairs_to_count_of_invs(pairs4);
		p4 = reversal(pairs_to_permut(pairs4));
		Permutation pp1(p1), pp2(p2);
		Permutation m = meet(pp1, pp2), jj = join(pp1, pp2);
		auto p5 = jj.get_value(), p6 = m.get_value();
		if (!weak_order_less(pp1, jj) || !weak_order_less(pp2, jj) || !weak_order_less(m, pp1) || !weak_order_less(m, pp2)) {
			cout << "Bad!" << endl;
		}
		if (p3 != p5 || p4 != p6) {
			cout << "Very Bad!" << endl;
		}
		if (i1 != pp1.cnt_invs || i2 != pp2.cnt_invs || i3 != jj.cnt_invs || i4 != m.cnt_invs) {
			cout << "Very Very Bad!" << endl;
		}
		print_vec(p5);
		print_vec(p1);
		print_vec(p2);
		print_vec(p6);
		cout << i1 << ' ' <<  i2 << ' ' << i3 << ' ' << i4 << endl;
		cout << pp1.cnt_invs << ' ' << pp2.cnt_invs << ' ' << jj.cnt_invs << ' ' << m.cnt_invs << endl;
	}
	*/

	short n = 4;
	int cnt = 0;
	vector<long double> parameters{0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999, 0.999999, 0.9999999, 1, 1.0000001, 1.000001, 1.00001, 1.0001, 1.001, 1.01, 1.02, 1.06, 1.1, 2};
	Diagram d = Diagram(n);
	vector<string> sn;
	for (auto it = d.graph->cbegin(); it != d.graph->cend(); ++it) {
		sn.push_back(it->first);
	}
	//reverse(sn.begin(), sn.end());
	vector<vector<string>> res;
	for (string p : sn) {
		TreeNode* node_p = d.graph->at(p);
		vector<vector<string>> tmp;
		tmp.push_back({ p });
		for (const vector<string>& v : res) {
			bool boo = true;
			for (const string& s : v) {
				TreeNode* ss = d.graph->at(s);
				if (weak_order_less(*(node_p->val), *(ss->val)) || weak_order_less(*(ss->val), *(node_p->val))) {
					boo = false;
					break;
				}
			}
			if (boo) {
				vector<string> to_add(v.cbegin(), v.cend());
				to_add.push_back(p);
				tmp.push_back(to_add);
			}
		}
		res.insert(res.end(), tmp.begin(), tmp.end());
		for (auto& v : tmp) {
			++cnt;
			vector<short> vec = d.increasing_set_inv_count(v);
			vector<long double> probabilities;
			for (long double q : parameters) {
				probabilities.push_back(mallows_measure_of_increasing_set(n, vec, q));
			}
			bool boo = true;
			for (int i = 1; i < probabilities.size(); ++i) {
				if (probabilities[i - 1] > probabilities[i]) {
					boo = false;
				}
			}
			if (!boo) {
				cout << endl;
				print_vec(probabilities);
				cout << endl;
			};
			cout << boo;
			if (cnt % 10000 == 0) cout << cnt/10000 << ' ';
		}
	}


	cout << cnt << endl;
	cout << "Bingo!" << endl;
	system("pause");
}
