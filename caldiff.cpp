#include "stdafx.h"
#include "iostream"

using namespace std;
double convert_to_julian(double, double, double, double, double, int);
void week(double);

int _tmain(int argc, _TCHAR* argv[])
{
	int  Seconds;
	double Year, Month, Day, Minutes, Hours, date1, date2, dif;

	cout << "Enter the first date (day month year hours minutes seconds):";
	cin >> Day
		>> Month
		>> Year
		>> Hours
		>> Minutes
		>> Seconds;

	date1 = convert_to_julian(Day, Month, Year, Hours, Minutes, Seconds);
	cout << date1 << endl;

	cout << "It's ";
	week(date1);

	cout << "Enter the second date (day month year hours minutes seconds):";
	cin >> Day
		>> Month
		>> Year
		>> Hours
		>> Minutes
		>> Seconds;

	date2 = convert_to_julian(Day, Month, Year, Hours, Minutes, Seconds);
	cout << date2 << endl;
	cout << "It's ";
	week(date2);

	dif = fabs(date1 - date2);

	cout << "Difference between these dates = " << dif << " days" << endl;

	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	cout.precision(5);

	system("pause");
	return 0;
}


double convert_to_julian(double Day, double Month, double Year, double Hours, double Minutes, int Seconds)
{
	double Mjd, MJD;
	int b, d;
	Minutes = Minutes + Seconds / 60;
	Hours = Hours + Minutes / 60;
	Day = Day + Hours / 24;
	if (Month <= 2) { Month += 12; --Year; }
	cout << "Choose Russia (enter 1) or Europe (enter 2): " << endl;
	cin >> d;
	if (d == 1)
	{
		if ((10000 * Year + 100 * Month + Day) <= 19183101)
			b = -2 + ((Year + 4716) / 4) - 1179;

		else b = (Year / 400) - (Year / 100) + (Year / 4);
		Mjd = 365 * Year - 679004 + b + int(30.6001*(Month + 1)) + Day;
	}
	else
		if ((10000 * Year + 100 * Month + Day) <= 15821004)
			b = -2 + ((Year + 4716) / 4) - 1179;
		else
			b = (Year / 400) - (Year / 100) + (Year / 4);
	Mjd = 365 * Year - 679004 + b + int(30.6001*(Month + 1)) + Day;
	MJD = 2400000.5 + Mjd;
	return MJD;
}

void week(double MJD)
{
	int w;
	w = int(MJD+0.5) % 7;

	switch (w) {
	case 0: {cout << "Monday" << endl;
		break; }
	case 1: {cout << "Tuesday" << endl;
		break; }
	case 2: {cout << "Wendesday" << endl;
		break; }
	case 3: {cout << "Thursday" << endl;
		break; }
	case 4: {cout << "Friday" << endl;
		break; }
	case 5: {cout << "Saturday" << endl;
		break; }
	case 6: {cout << "Sunday" << endl;
		break; }
	}
}

