#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#define EST_DIFF 1.0e-5

using namespace std;

void readMatrixFromFile(ifstream& input, int n, long double* a, long double* b, long double* c, long double* p, long double* q, long double* f, long double* r)
{
	//input >> n;

	//считывание массива p[] и f[0]
	for (int i = 0; i < n; i++)
		input >> p[i];
	input >> f[0];

	c[0] = p[n - 2];
	b[0] = p[n - 1];
	a[0] = 0;

	for (int i = 1; i < (n - 1); i++)
	{
		input >> c[i] >> b[i] >> a[i] >> f[i];
		r[i] = 0;
	}

	for (int i = 0; i < n; i++)
		input >> q[i];
	input >> f[n - 1];

	c[n - 1] = 0;
	b[n - 1] = q[0];
	a[n - 1] = q[1];


	//начальные значения массива r[]
	for (int i = 2; i < (n - 1); i++)
	{
		r[i] = 0;
	}
	r[0] = p[n - 1];
	r[1] = a[1];
	r[n - 1] = q[n - 1];
}

void generateMatrix(ofstream& file, int n, int deg)
{
	srand(time(NULL));
	long double rangeB = pow(10.0, deg);
	long double rangeA = (-1) * rangeB;

	file << n << endl;

	// 1-ая строка
	for (int i = 0; i < n + 1; ++i)
	{
		file << (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA << " ";
	}

	file << endl;

	// [1, n - 1] строки
	for (int i = 0; i < n - 2; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			file << (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA << " ";
		}

		file << endl;
	}

	// Последняя строка
	for (int i = 0; i < n + 1; ++i)
	{
		file << (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA << " ";
	}
}

void generateMatrixFromSingles(ofstream& file, int n, int deg, int mod1, int mod2)
{
	srand(time(NULL) + mod1 + mod2);
	long double rangeB = pow(10.0, deg);
	long double rangeA = (-1) * rangeB;
	long double summ = 0.0, arg = 0.0;

	file << n << endl;

	// 1-ая строка
	for (int i = 0; i < n; ++i)
	{
		arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
		file << arg << " ";
		summ += arg;
	}
	file << summ << " " << endl;
	summ = 0.0;

	// [1, n - 1] строки
	for (int i = 0; i < n - 2; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
			file << arg << " ";
			summ += arg;
		}
		file << summ << " " << endl;
		summ = 0.0;
	}

	// Последняя строка
	for (int i = 0; i < n; ++i)
	{
		arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
		file << arg << " ";
		summ += arg;
	}
	file << summ << " ";
}

void printMatrix(int n, long double* p, long double* q, long double* a, long double* b, long double* c, long double* f, long double* r)
{
	for (int i = 0; i < n; i++) //вывод первой строки
		cout << p[i] << "\t";
	cout << "|\t" << f[0] << endl;

	for (int i = 1; i < (n - 1); i++) //вывод строчек со второй по (n-1)-ую
	{
		for (int j = 0; j < (n - 3 - (i - 1)); j++)
			cout << "0\t";

		cout << c[i] << '\t' << b[i] << '\t' << a[i] << '\t';

		if (i > 1)
		{
			for (int j = 0; j < (i - 2); j++)
				cout << "0\t";
			cout << r[i] << '\t';
		}
		
		cout << "|\t" << f[i] << endl;
	}

	for (int i = 0; i < n; i++) //вывод последней строки
		cout << q[i] << '\t';
	cout << "|\t" << f[n - 1] << endl;
}

void printMatrixRounded(int n, long double* p, long double* q, long double* a, long double* b, long double* c, long double* f, long double* r)
{
	for (int i = 0; i < n; i++) //вывод первой строки
		cout << round(p[i] * 100) / 100 << "\t";
	cout << "|\t" << round(f[0] * 100) / 100 << endl;

	for (int i = 1; i < (n - 1); i++) //вывод строчек со второй по (n-1)-ую
	{
		for (int j = 0; j < (n - 3 - (i - 1)); j++)
			cout << "0\t";

		cout << round(c[i] * 100) / 100 << '\t' << round(b[i] * 100) / 100 << '\t' << round(a[i] * 100) / 100 << '\t';

		if (i > 1)
		{
			for (int j = 0; j < (i - 2); j++)
				cout << "0\t";
			cout << round(r[i] * 100) / 100 << '\t';
		}

		cout << "|\t" << round(f[i] * 100) / 100 << endl;
	}

	for (int i = 0; i < n; i++) //вывод последней строки
		cout << round(q[i] * 100) / 100 << '\t';
	cout << "|\t" << round(f[n - 1] * 100) / 100 << endl;
}

void solveSystem(int n, long double* p, long double* q, long double* a, long double* b, long double* c, long double* f, long double* r, long double* x)
{
	long double R;
	//STEP 1
	for (int i = 1; i < (n - 1); i++)
	{
		R = 1 / b[i];
		b[i] = 1;
		r[i] *= R;
		c[i] *= R;
		f[i] *= R;

		if (i == 1)
			a[1] = r[1];

		if (i != (n - 2))
		{
			R = a[i + 1];
			a[i + 1] = 0;
			r[i + 1] = (-1) * R * r[i];
			b[i + 1] -= R * c[i];
			f[i + 1] -= R * f[i];
		}

		R = p[n - i - 1];
		p[n - i - 1] = 0;
		p[n - 1] -= R * r[i];
		r[0] = p[n - 1];
		p[n - i - 2] -= R * c[i];
		f[0] -= R * f[i];

		R = q[n - i - 1];
		q[n - i - 1] = 0;
		q[n - 1] -= R * r[i];
		r[n - 1] = q[n - 1];
		q[n - i - 2] -= R * c[i];
		f[n - 1] -= R * f[i];
	}

	/*
	cout << endl;
	printMatrixRounded(n, p, q, a, b, c, f, r);
	cout << endl;
	*/

	//STEP 2
	R = 1 / p[n - 1];
	p[n - 1] = 1;
	r[0] = 1;
	p[0] *= R;
	f[0] *= R;

	R = q[n - 1];
	q[n - 1] = 0;
	r[n - 1] = 0;
	q[0] -= R * p[0];
	f[n - 1] -= R * f[0];

	R = 1 / q[0];
	q[0] = 1;
	f[n - 1] *= R;
	R = p[0];
	p[0] = 0;
	f[0] -= R * f[n - 1];

	/*
	cout << endl;
	printMatrixRounded(n, p, q, a, b, c, f, r);
	cout << endl;
	*/

	//STEP 3
	for (int i = 1; i < (n - 1); i++)
	{
		R = r[i];
		r[i] = 0;
		f[i] -= R * f[0];
	}
	a[1] = r[1];

	/*
	cout << endl;
	printMatrixRounded(n, p, q, a, b, c, f, r);
	cout << endl;
	*/

	//STEP 4 !!!
	x[n - 1] = f[0];
	x[0] = f[n - 1];
	for (int i = 1; i < (n - 1); i++)
	{
		x[i] = f[n - i - 1] - c[n - i - 1] * x[i - 1];
	}

	/*
	cout << "X: ";
	for (int i = 0; i < n; i++)
	{
		cout << x[i] << "  ";
	}
	cout << endl;
	*/
}

/*
bool checkResult(string fileName, int n, long double* a, long double* b, long double* c, long double* p, long double* q, long double* f, long double* r, long double* x)
{
	ifstream file(fileName);

	if (!file.is_open())
	{
		cerr << "Failed to open file " << fileName << endl;
		return -1;
	}

	file >> n;

	//long double* tempF = new long double[n];

	readMatrixFromFile(file, n, a, b, c, p, q, f, r);
	file.close();

	long double result = 0;

	// Проверка 0-ой строки
	for (int i = 0; i < n; ++i)
	{
		result += p[i] * x[i];
	}
	//cout << result << endl;

	if (result - f[0] > EST_DIFF)
	{
		return false;
	}

	result = 0;

	// Проверка n - 1 строки
	for (int i = 0; i < n; ++i)
	{
		result += q[i] * x[i];
	}
	//cout << result << endl;

	if (result - f[n - 1] > EST_DIFF)
	{
		return false;
	}

	// Проверка [1, n - 1] строк
	for (int i = 1; i < n - 1; ++i)
	{
		result = 0;

		result += c[i] * x[n - i - 2];
		result += b[i] * x[n - i - 1];
		result += a[i] * x[n - i];
		//cout << result << endl;

		if (result - f[i] > EST_DIFF)
		{
			return false;
		}
	}

	return true;
}
*/

/*
void findDiffFromSingles(int n, int m, long double* x)
{
	long double diff = 0;
	for (int i = 0; i < n; i++)
	{
		diff = max(diff, abs(x[i] - 1));
	}
	cout << n << 'x' << n << ", (-" << pow(10, m + 1) << ", " << pow(10, m + 1) << "):" << endl;
	cout << "Оценка точности решения = " << diff << endl;
}
*/

long double findDiffFromSingles(int n, int m, long double* x)
{
	long double diff = 0;
	for (int i = 0; i < n; i++)
	{
		diff = max(diff, abs(x[i] - 1));
	}
	
	return diff;
}

long double* generateRandMatrixRandX(ofstream& file, int n, int deg, int mod1, int mod2)
{
	srand(time(NULL) + mod1 + mod2);
	long double rangeB = pow(10.0, deg);
	long double rangeA = (-1) * rangeB;
	long double summ = 0.0, arg = 0.0;

	long double* genx = new long double[n];
	for (int i = 0; i < n; i++)
	{
		genx[i] = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
	}

	file << n << endl;

	// 1-ая строка
	for (int i = 0; i < n; ++i)
	{
		arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
		file << arg << " ";
		summ += arg * genx[i];
	}
	file << summ << " " << endl;
	summ = 0.0;

	// [1, n - 1] строки
	for (int i = 0; i < n - 2; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
			file << arg << " ";
			summ += arg * genx[n - i - 3 + j];
		}
		file << summ << " " << endl;
		summ = 0.0;
	}

	// Последняя строка
	for (int i = 0; i < n; ++i)
	{
		arg = (double)rand() / (double)RAND_MAX * (rangeB - rangeA) + rangeA;
		file << arg << " ";
		summ += arg * genx[i];
	}
	file << summ << " ";

	return genx;
}

/*
void findDiffSystem(int n, long double* acc_x, long double* x)
{
	long double diff = 0.0, curDiff = 0.0;
	for (int i = 0; i < n; i++)
	{
		curDiff = abs(x[i] - acc_x[i]);
		diff = max(diff, curDiff);
	}
	cout << "Средняя относительная погрешность системы = " << diff << endl;
}
*/

long double findDiffSystem(int n, int m, long double* acc_x, long double* x)
{
	long double diff = 0.0, curDiff = 0.0;
	for (int i = 0; i < n; i++)
	{
		curDiff = abs(x[i] - acc_x[i]);
		/*
		if (curDiff > pow(10, (m - 1)))
			curDiff /= acc_x[i];
		*/
		diff = max(diff, curDiff);
	}
	
	return diff;
}

void fullTest()
{
	string fileName;
	cout << "Введите имя файла для записей матрицы: ";
	cin >> fileName;
	cout << endl;

	long double differencesSingles[3][3] = { 0 };
	long double differencesSys[3][3] = { 0 };

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int numTest = 25;
			for (int k = 0; k < numTest; k++)
			{
				int n = i + 1;
				string testFile = fileName + to_string(i + 1) + to_string(j + 1) + ".txt";
				string testFileDiff = fileName + to_string(i + 1) + to_string(j + 1) + "Diff.txt";
				ofstream test(testFile);

				int degN = pow(10, n);
				generateMatrixFromSingles(test, degN, (j + 1), i, j);
				test.close();

				long double* p = new long double[degN];
				long double* q = new long double[degN];
				long double* f = new long double[degN];
				long double* a = new long double[degN];
				long double* b = new long double[degN];
				long double* c = new long double[degN];
				long double* r = new long double[degN];
				long double* x = new long double[degN];

				ifstream testRead(testFile);
				testRead >> degN;
				readMatrixFromFile(testRead, degN, a, b, c, p, q, f, r);
				solveSystem(degN, p, q, a, b, c, f, r, x);
				differencesSingles[i][j] += findDiffFromSingles(degN, j, x);

				ofstream testDiff(testFileDiff);
				long double* newx = new long double[i];
				newx = generateRandMatrixRandX(testDiff, degN, (j + 1), i, j);
				testDiff.close();

				ifstream testDiffRead(testFileDiff);
				testDiffRead >> degN;
				readMatrixFromFile(testDiffRead, degN, a, b, c, p, q, f, r);
				solveSystem(degN, p, q, a, b, c, f, r, x);
				differencesSys[i][j] += findDiffSystem(degN, j, newx, x);
			}
			//differencesSingles[i][j] /= numTest;
			//differencesSys[i][j] /= numTest;

			cout << pow(10, (i + 1)) << 'x' << pow(10, (i + 1)) << ", (-" << pow(10, j + 1) << ", " << pow(10, j + 1) << "):" << endl;
			cout << "Оценка точности решения = " << differencesSingles[i][j] / numTest << endl;
			cout << "Средняя относительная погрешность системы = " << differencesSys[i][j] / numTest << endl;
		}
		cout << endl;
	}
	
}

void accuracyAnswer(int n, long double* a, long double* b, long double* c, long double* p, long double* q, long double* r)
{
	long double* nf = new long double[n] {};
	long double* nx = new long double[n] {};



	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			for (int j = 0; j < n; j++)
			{
				nf[i] += p[j];
			}
		}
		else if (i == (n - 1))
		{
			for (int j = 0; j < n; j++)
			{
				nf[i] += q[j];
			}
		}
		else
		{
			nf[i] = a[i] + b[i] + c[i];
		}
	}

	
	solveSystem(n, p, q, a, b, c, nf, r, nx);

	long double diff = 0;
	for (int i = 0; i < n; i++)
	{
		diff = max(diff, abs(nx[i] - 1));
	}

	cout << "Оценка точности решения системы из файла равна " << diff << endl;

	delete[] nf;
	delete[] nx;
}

int main()
{
	setlocale(LC_ALL, "rus");
	
	string fileName;
	int n;

	cout << "1: Матрица из файла\t2: Сгенерировать тесты\nВыбор: ";
	cin >> n;

	if (n == 1)
	{
		cout << "Введите имя файла: ";
		cin >> fileName;
		fileName += ".txt";
		ifstream file(fileName);
		if (!file.is_open())
		{
			cout << "Не удалось открыть файл " << fileName << endl;
			return 0;
		}
		file >> n;

		long double* p = new long double[n];
		long double* q = new long double[n];
		long double* f = new long double[n];
		long double* a = new long double[n];
		long double* b = new long double[n];
		long double* c = new long double[n];
		long double* r = new long double[n];
		long double* x = new long double[n];

		readMatrixFromFile(file, n, a, b, c, p, q, f, r);

		solveSystem(n, p, q, a, b, c, f, r, x);

		cout << "X: ";
		for (int i = 0; i < n; i++)
		{
			cout << x[i] << "  ";
		}
		cout << endl;

		readMatrixFromFile(file, n, a, b, c, p, q, f, r);
		accuracyAnswer(n, a, b, c, p, q, r);

		delete[] p;
		delete[] q;
		delete[] f;
		delete[] a;
		delete[] b;
		delete[] c;
		delete[] r;
		delete[] x;

		return 0;
	}
	
	else if (n == 2)
	{
		fullTest();
	}

	else
	{
		cout << "Неверный номер ответа.";
		return 0;
	}
}