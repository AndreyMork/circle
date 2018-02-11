#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <map>


using namespace std;


struct Point
{
	double x;
	double y;
	int index; // required for deleting
};


void printPoints(vector <Point*> &points)
{
	int c = 1;
	for (Point* i : points)
	{
		cout << c++ << ". (" << i->x << ", " << i->y << ")" << endl;
	}
}


void generatePoints(int n, vector <Point*> &points)
{
	uniform_real_distribution<double> dist(-1000, 1000);  //(min, max)
	mt19937 rng;
	rng.seed(std::random_device{}());

	for (int i = 0; i < n; i++)
	{
		(points).push_back(new Point{ dist(rng), dist(rng) });
	}
}


void getRminmax(vector <Point*> &points, double& rmin, double& rmax)
{
	double Ymax = abs(points[0]->y);
	rmin = Ymax;
	double tempY;

	double xl = points[0]->x;
	double xr = xl;

	for (Point* i : points)
	{
		tempY = abs(i->y);
		if (tempY > Ymax)
			Ymax = tempY;
		else if (tempY < rmin)
			rmin = tempY;

		if (i->x < xl)
			xl = i->x;
		else if (i->x > xr)
			xr = i->x;
	}

	rmax = sqrt((xr - xl) * (xr - xl) + Ymax * Ymax);
}


struct Bound
{
	double x;
	bool begin;
	Point* origin;

	bool operator()(const Bound* b1, const Bound* b2) const
	{
		if (b1->x == b2->x)
		{
			if (b1->origin == b2->origin) {
				return b1->begin > b2->begin;
			}
			else
			{
				return b1->begin < b2->begin;
			}
		}
		else
		{
			return b1->x < b2->x;
		}
	}
};

void printBounds(vector <Bound*> &bounds)
{
	int counter = 0;
	for (Bound* b : bounds)
	{
		counter += (b->begin ? 1 : -1);
		cout << counter << "   " << b->x << "   " << b->begin << "   " << b->origin << endl;
	}
}


void getBounds(vector <Point*> &points, double r, vector <Bound*> &bounds)
{
	bounds.clear();
	double temp;

	for (Point* p : points)
	{
		temp = r * r - p->y * p->y;
		if (temp >= 0)
		{
			bounds.push_back(new Bound{ p->x - sqrt(temp), true, p });
			bounds.push_back(new Bound{ p->x + sqrt(temp), false, p });
		}
	}

	sort(bounds.begin(), bounds.end(), Bound()); // O(nlogn)
}


void getNewPoints(vector <Bound*> &bounds, vector<Point*> &newpoints, int k)
{
	newpoints.clear();
	if (bounds.size() < 2 * k) 	{return;}

	int counter = 0, k_counter = 0;
	
	int index = 0;
	bool k_reached = false;
	for (Bound* b : bounds)
	{
		if (b->begin)
		{
			counter++;
			if (counter >= k)
			{
				k_reached = true;
				k_counter++;
			}
		}
		else
		{
			counter--;
			if (k_reached)
			{
				newpoints.push_back(b->origin);
				b->origin->index = index++;
			}
			if (counter == 0)
			{
				k_reached = false;
			}
		}
	}

	if (k_counter == 1)
	{
		counter = 0;
		k_reached = false;
		for (Bound* b : bounds)
		{
			if (b->begin)
			{
				if (!k_reached)
					{
						counter++;
						if (counter == k)
						{
							k_reached = true;
						}
					}
				else
				{
					newpoints[b->origin->index] = nullptr; 	// delete
				}
			}
			else
			{
				counter--;
				if (counter == 0)
				{
					k_reached = false;
				}
			}
		}
		newpoints.erase(remove(newpoints.begin(), newpoints.end(), nullptr), newpoints.end());
	}
}



#define N 1000
#define K 2
#define Delta 10e-6
int main()
{
	vector <Point*> inpoints;
	generatePoints(N, inpoints);

	vector <Point*> points = inpoints; // to not mess up the original data

	points.push_back(new Point{ 0, 1 });
	points.push_back(new Point{ 2, 1 });

	double rmax, rmin, r = 0, rl = 1;	
	getRminmax(points, rmin, rmax);

	vector<Bound*> bounds;
	vector<Point*> newpoints;

	int counter = 0;
	while (abs(rl - r) > Delta || points.size() < K)
	{
		rl = r;
		r = (rmax + rmin) / 2;

		getBounds(points, r, bounds);
		getNewPoints(bounds, newpoints, K);

		if (newpoints.size() < K)
		{
			rmin = r;
		}
		else if (newpoints.size() >= K)
		{
			rmax = r;
			points = newpoints;
		}

		if (bounds.size() < 2 * K)
			rl = r + 1;	// required to get centre right

		// cout << abs(r - rl) << endl;
		cout << ++counter << "." << endl;
		cout << "p = " << newpoints.size() << endl;
		cout << "r = " << r << endl << endl;
		// printBounds(bounds);
		// cin.get();
	}

	cout << endl << "-----------" << endl;
	getBounds(points, r, bounds);
	// printBounds(bounds);
	// cout << endl;
	// cout << bounds.size() << endl;
	cout << "c in [" << bounds[K - 1]->x << ", " << bounds[K]->x << "]" << endl;
	double centre = (bounds[K]->x + bounds[K - 1]->x) / 2;


	cout << "r = " << r << endl;
	cout << "c = " << centre << endl;
	int c = 1;
	for (Point* p : points)
	{
		cout << c++ << ". " << "(" << p->x << ", " << p->y << ")\n    " << (p->x - centre) * (p->x - centre) + p->y * p->y << " <= " << (r + Delta) * (r + Delta)
			<< "\t" << ((p->x - centre) * (p->x - centre) + p->y * p->y <= (r + Delta) * (r + Delta) ? "true" : "false") << endl;
	}
	
	cout << "...";
	cin.get();
	return 0;
}