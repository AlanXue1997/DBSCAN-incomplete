#include <cstdio>
#include <time.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Windows.h>
#include <stack>
#include <algorithm>

//#define DEBUG
//#define MULTI	//Won't process missing data

#define pointSet std::vector<Point*>

#define MAX_FIND 100
#define epsilon 0.0000001
#define sample_n 1000
#define NOT_VISITED -1
#define MAX_CLUSTER 100
#define INFINITE 10000

#define ATTRIBUTIONS 3
#define CLUSTER 2	//Number of original clusters / Number of label types

#define K 4

class Point {
public:
	double a[ATTRIBUTIONS];
	Point() {
		ddrs = NULL;
		flag = NOT_VISITED;
		incomplete = 0;
	}
	int flag;
	int flag_ans;
	int incomplete;
	pointSet *ddrs;
};

class TREE_Point {
public:
	Point *point;
	int d;
	TREE_Point *l, *r;
};

double distance(Point *a, Point *b) {
	double ans = 0;
	for (int i = 0; i < ATTRIBUTIONS; ++i) {
		ans += (a->a[i] - b->a[i])*(a->a[i] - b->a[i]);
	}
	return sqrt(ans);
}

double pseDistance(Point *a, Point *b) {
	double ans = 0;
	for (int i = 0; i < ATTRIBUTIONS; ++i) {
		if(b->a[i] != INFINITE) ans += (a->a[i] - b->a[i])*(a->a[i] - b->a[i]);
	}
	return sqrt(ans);
}

class KD_TREE {
private:
	TREE_Point* root;
	double E;
	int findSplit(pointSet *set) {
		double var[ATTRIBUTIONS];
		double sum[ATTRIBUTIONS];
		memset(var, 0, sizeof(var));
		memset(sum, 0, sizeof(sum));
		for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
			for (int j = 0; j < ATTRIBUTIONS; ++j) {
				sum[j] += (*i)->a[j];
			}
		}
		for (int j = 0; j < ATTRIBUTIONS; ++j) {
			sum[j] /= set->size();
		}
		for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
			for (int j = 0; j < ATTRIBUTIONS; ++j) {
				var[j] += ((*i)->a[j] - sum[j])*((*i)->a[j] - sum[j]);
			}
		}
		int k = 0;
		for (int i = 1; i < ATTRIBUTIONS; ++i) {
			if (var[k] < var[i]) k = i;
		}
		for (int i = 1; i < set->size(); ++i) {
			int j = i;
			while (j > 0 && set->at(j - 1)->a[k] > set->at(j)->a[k]) {
				std::swap(set->at(j - 1), set->at(j));
				j--;
			}
		}
		return k;
	}

	void findNeighbor(pointSet *set, TREE_Point *p, Point *point) {
		double dis = distance(p->point, point);
		if (dis < E && dis > epsilon) set->push_back(p->point);
		if (point->a[p->d] < p->point->a[p->d]) {
			if (p->l != NULL) findNeighbor(set, p->l, point);
			if (p->r != NULL && (p->point->a[p->d] - point->a[p->d]) < E) findNeighbor(set, p->r, point);
		}
		else {
			if (p->r != NULL) findNeighbor(set, p->r, point);
			if (p->l != NULL && (point->a[p->d] - p->point->a[p->d]) < E) findNeighbor(set, p->l, point);
		}
	}

	TREE_Point* construct(pointSet *set) {
		TREE_Point* p = new TREE_Point();
		if (set->size() == 0) return NULL;
		if (set->size() == 1) {
			p->d = -1;
			p->l = p->r = NULL;
			p->point = set->at(0);
			return p;
		}
		pointSet* sample = new pointSet();
		if (set->size() <= sample_n) {
			sample = set;
		}
		else {
			for (int i = 0; i < set->size(); ++i) {
				if (rand() % (set->size() - i) < (sample_n - sample->size())) {
					sample->push_back(set->at(i));
				}
				if (sample->size() == sample_n) break;
			}
		}
		p->d = findSplit(sample);
		p->point = sample->at(sample->size() / 2);
		pointSet* lset;
		pointSet* rset;
		if (set->size() <= sample_n)
		{
			lset = new pointSet(set->size() / 2);
			std::copy(set->begin(), set->begin() + set->size() / 2, lset->begin());
			if (set->size() > 2) {
				rset = new pointSet(set->size() - set->size() / 2 - 1);
				std::copy(set->begin() + set->size() / 2 + 1, set->end(), rset->begin());
			}
			else {
				rset = NULL;
			}
			delete set;
			p->l = construct(lset);
			if (rset != NULL) {
				p->r = construct(rset);
			}
			else {
				p->r = NULL;
			}
		}
		else {
			delete sample;
			lset = new pointSet();
			rset = new pointSet();
			int mark = 0;
			for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
				int &k = p->d;
				if ((*i)->a[k] < p->point->a[k]) {
					lset->push_back(*i);
				}
				else if ((*i)->a[k] > p->point->a[k]) {
					rset->push_back(*i);
				}
				else {
					if (mark) {
						lset->push_back(*i);
					}
					else {
						rset->push_back(*i);
					}
					mark = !mark;
				}
			}
			delete set;
			p->l = construct(lset);
			p->r = construct(rset);
		}
		return p;
	}

	void findPredict(std::vector<double> *left, std::vector<double> *right, TREE_Point *p, Point *point, int d) {
		if (p == NULL) return;
		if (left->size() > MAX_FIND) return;
		double dis = pseDistance(p->point, point);
		if (dis < E) {
			dis = sqrt(E*E - dis*dis);
			left->push_back(p->point->a[d] - dis); std::push_heap(left->begin(), left->end());
			right->push_back(p->point->a[d] + dis); std::push_heap(right->begin(), right->end());
		}
		if (point->a[p->d] != INFINITE) {
			if (point->a[p->d] < p->point->a[p->d]) {
				findPredict(left, right, p->l, point, d);
				if(p->point->a[p->d] - point->a[p->d] < E) findPredict(left, right, p->r, point, d);
			}
			else {
				findPredict(left, right, p->r, point, d);
				if (point->a[p->d] - p->point->a[p->d] < E) findPredict(left, right, p->l, point, d);
			}
		}
		else {
			findPredict(left, right, p->r, point, d);
			findPredict(left, right, p->l, point, d);
		}
	}

public:
	KD_TREE(pointSet *set) {
		pointSet *set_copy = new pointSet(set->size());
		std::copy(set->begin(), set->end(), set_copy->begin());
		root = new TREE_Point();
		root->r = construct(set_copy);
	}

	pointSet* neighbor(Point *point, double E) {
		pointSet *set = new pointSet();
		this->E = E;
		if (root->r != NULL) findNeighbor(set, root->r, point);
		return set;
	}

	double predict(Point *point, double E, int d) {
		this->E = E;
		std::vector<double> *left = new std::vector<double>();
		std::vector<double> *right = new std::vector<double>();
		findPredict(left, right, root->r, point, d);
		int k = 0, max = 0;
		double l, r;
		while (!left->empty()) {
			if (!right->empty() && left->front() <= right->front()) {
				k++;
				if (k >= max) {
					r = right->front();
					max = k;
				}
				std::pop_heap(right->begin(), right->end()); right->pop_back();
			}
			else {
				if (k == max) l = left->front();
				k--;
				std::pop_heap(left->begin(), left->end()); left->pop_back();
			}
		}
		return (l + r) / 2;
	}

	void insert(Point *point, double E) {
		TREE_Point *p = root->r;
		if (p == NULL) {
			root->r = new TREE_Point();
			root->r->l = root->r->r = NULL;
			root->d = 0;
			root->point = point;
			return;
		}
		while (1) {
			if (point->a[p->d] < p->point->a[p->d]) {
				if (p->l == NULL) {
					p->l = new TREE_Point();
					p = p->l;
					p->d = rand() % ATTRIBUTIONS;
					p->l = p->r = NULL;
					p->point = point;
					return;
				}
				else {
					p = p->l;
				}
			}
			else {
				if (p->r == NULL) {
					p->r = new TREE_Point();
					p = p->r;
					p->d = rand() % ATTRIBUTIONS;
					p->l = p->r = NULL;
					p->point = point;
					return;
				}
				else {
					p = p->r;
				}
			}
		}
	}

#ifdef DEBUG
	void _traverse(TREE_Point *p) {
		if (p == NULL) {
			std::cout << "*" << std::endl;
			return;
		}
		for (int i = 0; i < ATTRIBUTIONS; ++i) {
			std::cout << p->point->a[i] << ' ';
		}
		std::cout << std::endl;
		_traverse(p->l);
		_traverse(p->r);
	}
	void traverse() {
		if (root->r == NULL) std::cout << "(Empty)" << std::endl;
		else _traverse(root->r);
	}
#endif
};

void DFS(pointSet *set, Point* x, int flag) {
	std::stack<Point*> q;
	int k = 1;
	q.push(x);
	x->flag = flag;
	while (!q.empty()) {
		Point* p = q.top();
		q.pop();
		for (pointSet::iterator i = p->ddrs->begin(); i != p->ddrs->end(); ++i) {
			if ((*i)->flag == NOT_VISITED) {
				(*i)->flag = flag;
				k++;
				if ((*i)->ddrs != NULL) q.push(*i);
			}
		}
	}
}

void process(pointSet *set, double E) {
	int flag = 0;
#ifdef MULTI
	for (pointSet::iterator i = set->begin(); i != set->end(); i++) {
		(*i)->flag = NOT_VISITED;
	}
#endif
#ifndef MULTI
	std::cout << "Clustering...";
#endif
	for (int i = 0; i < set->size(); ++i) {
		if (set->at(i)->flag == NOT_VISITED && set->at(i)->ddrs != NULL) {
			DFS(set, set->at(i), flag);
			flag++;
		}
	}
#ifndef MULTI
	std::cout << "[Done]" << std::endl;
#endif
}

void normalize(pointSet *set) {
	double s_max[ATTRIBUTIONS];
	double s_min[ATTRIBUTIONS];
	double num[ATTRIBUTIONS] = { 0 };
	double avg[ATTRIBUTIONS] = { 0 };
	for (int i = 0; i < ATTRIBUTIONS; ++i) {
		s_max[i] = -INFINITE;
		s_min[i] = INFINITE;
	}
	for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
		for (int j = 0; j < ATTRIBUTIONS; ++j) {
			if ((*i)->a[j] != INFINITE) {
				num[j]++;
				if (s_max[j] < (*i)->a[j]) s_max[j] = (*i)->a[j];
				if (s_min[j] > (*i)->a[j]) s_min[j] = (*i)->a[j];
				avg[j] += (*i)->a[j];
			}
		}
	}
	for (int i = 0; i < ATTRIBUTIONS; ++i) avg[i] /= num[i];
	for (pointSet::iterator p = set->begin(); p != set->end(); ++p) {
		for (int i = 0; i < ATTRIBUTIONS; ++i) {
			if ((*p)->a[i] != INFINITE) {
				(*p)->a[i] = ((*p)->a[i] - avg[i]) / (s_max[i] - s_min[i]);
			}
		}
	}
}

double accuracy(pointSet *set) {
	int cluster[MAX_CLUSTER][CLUSTER];
	int a_num[CLUSTER];
	int cluster_num[MAX_CLUSTER];
	int flag_max = 0;
	memset(cluster, 0, sizeof(cluster));
	memset(a_num, 0, sizeof(a_num));
	memset(cluster_num, 0, sizeof(cluster_num));

	for (pointSet::iterator i = set->begin(); i != set->end(); i++) {
		if ((*i)->flag != NOT_VISITED && (*i)->flag_ans != NOT_VISITED) {
			cluster[(*i)->flag][(*i)->flag_ans]++;
			a_num[(*i)->flag_ans]++;
			cluster_num[(*i)->flag]++;
			if (flag_max < (*i)->flag) flag_max = (*i)->flag;
		}
	}

	double ans = 0;
	for (int i = 0; i <= flag_max; ++i) {
		double fi = 0;
		for (int j = 0; j < CLUSTER; ++j) {
			double p = (double)cluster[i][j] / (double)cluster_num[i];
			double r = (double)cluster[i][j] / (double)a_num[j];
			double f = 2 * p*r / (p + r);
			if (fi < f) fi = f;
		}
		ans += ((double)cluster_num[i] / (double)set->size())*fi;
	}
	return ans;
}

int main(int argc, char* argv[])
{

	srand(1);
	double E;
	int N;
	double l, r, m;
	std::fstream fin;
#ifdef MULTI
	std::cin >> l >> r;
	fin.open("input.txt", std::fstream::in);
#else
	if (argc > 1) {
		E = atof(argv[2]);
		fin.open(argv[1], std::fstream::in);
	}
	else
	{
		fin.open("input.txt", std::fstream::in);
		std::cin >> E;
	}
	LARGE_INTEGER nFreq, nNow, nRead, nDone;
	double t_read, t_cal, t_all;
	QueryPerformanceFrequency(&nFreq);
	QueryPerformanceCounter(&nNow);
#endif
	pointSet* set = new pointSet();
	pointSet* completeSet;// = new pointSet();
	pointSet* incompleteSet[ATTRIBUTIONS+1];
	pointSet* missingSet = new pointSet();
	for (int i = 0; i <= ATTRIBUTIONS; ++i) {
		incompleteSet[i] = new pointSet();
	}
	//std::ifstream fin("input.txt");
	std::cout << "Scanning...";
	fin >> N;
	for (int i = 0; i < N; ++i) {
		char c;
		Point* p = new Point();
		for (int j = 0; j < ATTRIBUTIONS; ++j) {
			if (!(fin >> p->a[j])) {
				fin.clear();
				do { fin >> c; } while (c != '?');
				p->incomplete++;
				p->a[j] = INFINITE;
			}
		}
		incompleteSet[p->incomplete]->push_back(p);
		if (!(fin >> p->flag_ans)) {
			fin.clear();
			do { fin >> c; } while (c != '?');
			p->flag_ans = NOT_VISITED;
		}
		/*
		if (p->incomplete == 0) {
			completeSet->push_back(p);
		}
		else {
			missingSet->push_back(p);
		}*/
		set->push_back(p);
	}
	completeSet = incompleteSet[0];
	fin.close();
	std::cout << "[Done]" << std::endl;
#ifndef MULTI
	QueryPerformanceCounter(&nRead);
#endif
	std::cout << "Normalizing...";
	normalize(set);
	std::cout << "[Done]" << std::endl;
	std::cout << "Constructing KD-Tree...";
	KD_TREE tree(completeSet);
	std::cout << "[Done]" << std::endl;
#ifdef MULTI
	for (double E = l; E < r; E += (r - l) / 10) {
		for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
			delete (*i)->ddrs;
			(*i)->ddrs = tree.neighbor(*i, E);
			if ((*i)->ddrs->size() < K) (*i)->ddrs = NULL;
		}
		process(set, E);
		std::cout << "E = " << E << ", accuracy = " << accuracy(set) << std::endl;
	}
#else
	/*
	int d;
	do {
		d = ATTRIBUTIONS;
		for (pointSet::iterator p = missingSet->begin(); p != missingSet->end(); p++) {
			if ((*p)->incomplete > 0 && (*p)->incomplete < d) d = (*p)->incomplete;
		}
		for (pointSet::iterator p = missingSet->begin(); p != missingSet->end(); p++) {
			if ((*p)->incomplete == d) {
				for (int i = 0; i < ATTRIBUTIONS; ++i) {
					if ((*p)->a[i] == INFINITE) {
						(*p)->a[i] = tree.predict((*p), E, i);
						(*p)->incomplete--;
						if ((*p)->incomplete == 0) tree.insert(*p, E);
					}
				}
			}
		}
	} while (d < ATTRIBUTIONS);*/
	int sum = 0;
	for (int i = 1; i <= ATTRIBUTIONS; ++i) sum += incompleteSet[i]->size();
	std::cout << "Processing incomplete data: " << sum << "in total" << std::endl;
	for (int k = 1; k <= ATTRIBUTIONS; ++k) {
		for (pointSet::iterator p = incompleteSet[k]->begin(); p != incompleteSet[k]->end(); p++) {
			for (int i = 0; i < ATTRIBUTIONS; ++i) {
				if ((*p)->a[i] == INFINITE) {
					(*p)->a[i] = tree.predict((*p), E, i);
				}
			}
			tree.insert(*p, E);
#ifdef DEBUG
			std::cout << sum-- << std::endl;
#endif
		}
		delete incompleteSet[k];
	}
	
	std::cout << "Range searching...";
	for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
		(*i)->ddrs = tree.neighbor(*i, E);
		if ((*i)->ddrs->size() < K) {
			delete (*i)->ddrs;
			(*i)->ddrs = NULL;
		}
	}
	std::cout << "[Done]" << std::endl;
	process(set, E);

	std::ofstream fout("output.txt");
	for (pointSet::iterator i = set->begin(); i != set->end(); ++i) {
		fout << (*i)->flag_ans << '\t' << (*i)->flag << std::endl;
	}
	fout.close();

	double acc = accuracy(set);
	std::cout << "accuracy = " << acc << std::endl;
	QueryPerformanceCounter(&nDone);
	t_read = (double)(nRead.QuadPart - nNow.QuadPart) / (double)nFreq.QuadPart;
	t_cal = (double)(nDone.QuadPart - nRead.QuadPart) / (double)nFreq.QuadPart;
	t_all = t_read + t_cal;
	std::cout << "Reading time:\t" << t_read << std::endl;
	std::cout << "Calculate time:\t" << t_cal << std::endl;
	std::cout << "Total time:\t" << t_all << std::endl;
	std::ofstream log_out("log.txt", std::ios::app);
	log_out << acc << ' ' << t_cal << std::endl;
#endif
	if (argc <= 1) {
		system("pause");
	}
	return 0;
}
