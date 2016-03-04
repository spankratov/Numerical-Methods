#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

#include "armadillo"

using namespace std;
using namespace arma;

typedef double(*funn)(const vec &arg);

double firstTask1(const vec& arg) { return cos(arg[0] + 0.5) + arg[1] - 0.8; }
double firstTask2(const vec& arg) { return sin(arg[1]) - 2 * arg[0] - 1.6; }
double firstTaskJacoby11(const vec& arg) { return -sin(arg[0] + 0.5); }
double firstTaskJacoby12(const vec& arg) { return 1; }
double firstTaskJacoby21(const vec& arg) { return -2; }
double firstTaskJacoby22(const vec& arg) { return cos(arg[1]); }
vector<funn> firstTaskVec = { firstTask1, firstTask2 };
vector< vector<funn> > firstTaskJacoby = { { firstTaskJacoby11, firstTaskJacoby12 }, { firstTaskJacoby21, firstTaskJacoby22 } };

/*double secondTask1(const vec& arg) { return cos(arg[0] * arg[1]) - exp(-3 * arg[2]) + arg[3] * arg[4] * arg[4] - arg[5] - sinh(2 * arg[7]) * arg[8] + 2 * arg[9] + 2.0004339741653854440; }
double secondTask2(const vec& arg) { return sin(arg[0] * arg[1]) + arg[2] * arg[8] * arg[6] - exp(-arg[9] + arg[5]) + 3 * arg[4] * arg[4] - arg[5] * (arg[7] + 1) + 10.886272036407019994; }
double secondTask3(const vec& arg) { return arg[0] - arg[1] + arg[2] - arg[3] + arg[4] - arg[5] + arg[6] - arg[7] + arg[8] - arg[9] - 3.1361904761904761904; }
double secondTask4(const vec& arg) { return 2 * cos(-arg[8] + arg[3]) + arg[4] / (arg[2] + arg[0]) - sin(arg[1] * arg[1]) + pow(cos(arg[6] * arg[9]), 2) - arg[7] - 0.1707472705022304757; }
double secondTask5(const vec& arg) { return sin(arg[4]) + 2 * arg[7] * (arg[2] + arg[0]) - exp(-arg[6] * (-arg[9] + arg[5])) + 2 * cos(arg[1]) - 1 / (arg[3] - arg[8]) - 0.3685896273101277862; }
double secondTask6(const vec& arg) { return exp(arg[0] - arg[3] - arg[8]) + arg[4] * arg[4] / arg[7] + 0.5*cos(3 * arg[9] * arg[1]) - arg[5] * arg[2] + 2.0491086016771875115; }
double secondTask7(const vec& arg) { return arg[1] * arg[1] * arg[1] * arg[6] - sin(arg[9] / arg[4] + arg[7]) + (arg[0] - arg[5])*cos(arg[3]) + arg[2] - 0.7380430076202798014; }
double secondTask8(const vec& arg) { return arg[4] * pow(arg[0] - 2 * arg[5], 2) - 2 * sin(-arg[8] + arg[2]) + 1.5*arg[3] - exp(arg[1] * arg[6] + arg[9]) + 3.5668321989693809040; }
double secondTask9(const vec& arg) { return 7 / arg[5] + exp(arg[4] + arg[3]) - 2 * arg[1] * arg[7] * arg[9] * arg[6] + 3 * arg[8] - 3 * arg[0] - 8.4394734508383257499; }
double secondTask10(const vec& arg) { return arg[9] * arg[0] + arg[8] * arg[1] - arg[7] * arg[2] + sin(arg[3] + arg[4] + arg[5])*arg[6] - 0.78238095238095238096; }
double secondTaskJacoby11(const vec& arg) { return -sin(arg[0] * arg[1])*arg[1]; }
double secondTaskJacoby12(const vec& arg) { return -sin(arg[0] * arg[1])*arg[0]; }
double secondTaskJacoby13(const vec& arg) { return 3 * exp(-3 * arg[2]); }
double secondTaskJacoby14(const vec& arg) { return arg[4] * arg[4]; }
double secondTaskJacoby15(const vec& arg) { return 2 * arg[3] * arg[4]; }
double secondTaskJacoby16(const vec& arg) { return -1; }
double secondTaskJacoby17(const vec& arg) { return 0; }
double secondTaskJacoby18(const vec& arg) { return -2 * cosh(2 * arg[7])*arg[8]; }
double secondTaskJacoby19(const vec& arg) { return -sinh(2 * arg[7]); }
double secondTaskJacoby110(const vec& arg) { return 2; }
double secondTaskJacoby21(const vec& arg) { return cos(arg[0] * arg[1])*arg[1]; }
double secondTaskJacoby22(const vec& arg) { return cos(arg[0] * arg[1])*arg[0]; }
double secondTaskJacoby23(const vec& arg) { return arg[8] * arg[6]; }
double secondTaskJacoby24(const vec& arg) { return 0; }
double secondTaskJacoby25(const vec& arg) { return 6 * arg[4]; }
double secondTaskJacoby26(const vec& arg) { return -exp(-arg[9] + arg[5]) - arg[7] - 1; }
double secondTaskJacoby27(const vec& arg) { return arg[2] * arg[8]; }
double secondTaskJacoby28(const vec& arg) { return -arg[5]; }
double secondTaskJacoby29(const vec& arg) { return arg[2] * arg[6]; }
double secondTaskJacoby210(const vec& arg) { return exp(-arg[9] + arg[5]); }
double secondTaskJacoby31(const vec& arg) { return 1; }
double secondTaskJacoby32(const vec& arg) { return -1; }
double secondTaskJacoby33(const vec& arg) { return 1; }
double secondTaskJacoby34(const vec& arg) { return -1; }
double secondTaskJacoby35(const vec& arg) { return 1; }
double secondTaskJacoby36(const vec& arg) { return -1; }
double secondTaskJacoby37(const vec& arg) { return 1; }
double secondTaskJacoby38(const vec& arg) { return -1; }
double secondTaskJacoby39(const vec& arg) { return 1; }
double secondTaskJacoby310(const vec& arg) { return -1; }
double secondTaskJacoby41(const vec& arg) { return -arg[4] / pow(arg[2] + arg[0], 2); }
double secondTaskJacoby42(const vec& arg) { return -2 * cos(arg[1] * arg[1])*arg[1]; }
double secondTaskJacoby43(const vec& arg) { return -arg[4] / pow(arg[2] + arg[0], 2); }
double secondTaskJacoby44(const vec& arg) { return -2 * sin(-arg[8] + arg[3]); }
double secondTaskJacoby45(const vec& arg) { return 1 / (arg[2] + arg[0]); }
double secondTaskJacoby46(const vec& arg) { return 0; }
double secondTaskJacoby47(const vec& arg) { return -2 * cos(arg[6] * arg[9])*sin(arg[6] * arg[9])*arg[9]; }
double secondTaskJacoby48(const vec& arg) { return -1; }
double secondTaskJacoby49(const vec& arg) { return 2 * sin(-arg[8] + arg[3]); }
double secondTaskJacoby410(const vec& arg) { return -2 * cos(arg[6] * arg[9])*sin(arg[6] * arg[9])*arg[6]; }
double secondTaskJacoby51(const vec& arg) { return 2 * arg[7]; }
double secondTaskJacoby52(const vec& arg) { return -2 * sin(arg[1]); }
double secondTaskJacoby53(const vec& arg) { return 2 * arg[7]; }
double secondTaskJacoby54(const vec& arg) { return pow(-arg[8] + arg[3], -2); }
double secondTaskJacoby55(const vec& arg) { return cos(arg[4]); }
double secondTaskJacoby56(const vec& arg) { return arg[6] * exp(-arg[6] * (-arg[9] + arg[5])); }
double secondTaskJacoby57(const vec& arg) { return -(arg[9] - arg[5])*exp(-arg[6] * (-arg[9] + arg[5])); }
double secondTaskJacoby58(const vec& arg) { return 2 * arg[2] + arg[0]; }
double secondTaskJacoby59(const vec& arg) { return -1 / pow(-arg[8] + arg[3], 2); }
double secondTaskJacoby510(const vec& arg) { return -arg[6] * exp(-arg[6] * (-arg[9] + arg[5])); }

double secondTaskJacoby61(const vec& arg) {
	return exp(arg[0] - arg[3] - arg[8]);
}
double secondTaskJacoby62(const vec& arg) {
	return -3 / 2 * sin(3 * arg[9] * arg[1])*arg[9];
}
double secondTaskJacoby63(const vec& arg) {
	return -arg[5];
}
double secondTaskJacoby64(const vec& arg) {
	return -exp(arg[0] - arg[3] - arg[8]);
}
double secondTaskJacoby65(const vec& arg) {
	return 2 * arg[4] / arg[7];
}
double secondTaskJacoby66(const vec& arg) {
	return -arg[2];
}
double secondTaskJacoby67(const vec& arg) {
	return 0;
}
double secondTaskJacoby68(const vec& arg) {
	return -(arg[4] * arg[4]) / (arg[7] * arg[7]);
}
double secondTaskJacoby69(const vec& arg) {
	return -exp(arg[0] - arg[3] - arg[8]);
}
double secondTaskJacoby610(const vec& arg) {
	return -3 / 2 * sin(3 * arg[9] * arg[1]) * arg[1];
}
double secondTaskJacoby71(const vec& arg) {
	return cos(arg[3]);
}
double secondTaskJacoby72(const vec& arg) {
	return 3 * arg[1] * arg[1] * arg[6];
}
double secondTaskJacoby73(const vec& arg) {
	return 1;
}
double secondTaskJacoby74(const vec& arg) {
	return -(arg[0] - arg[5])*sin(arg[3]);
}
double secondTaskJacoby75(const vec& arg) {
	return cos(arg[9] / arg[4] + arg[7])*arg[9] * pow(arg[4], -2);
}
double secondTaskJacoby76(const vec& arg) {
	return -cos(arg[3]);
}
double secondTaskJacoby77(const vec& arg) {
	return arg[1] * arg[1];
}
double secondTaskJacoby78(const vec& arg) {
	return -cos(arg[9] / arg[4] + arg[7]);
}
double secondTaskJacoby79(const vec& arg) {
	return 0;
}
double secondTaskJacoby710(const vec& arg) {
	return -cos(arg[9] / arg[4] + arg[7]) / arg[4];
}
double secondTaskJacoby81(const vec& arg) {
	return 2 * arg[4] * (arg[0] - 2 * arg[5]);
}
double secondTaskJacoby82(const vec& arg) {
	return -arg[6] * exp(arg[1] * arg[6] + arg[9]);
}
double secondTaskJacoby83(const vec& arg) {
	return -2 * cos(-arg[8] + arg[2]);
}
double secondTaskJacoby84(const vec& arg) {
	return 1.5;
}
double secondTaskJacoby85(const vec& arg) {
	return pow(arg[0] - 2 * arg[5], 2);
}
double secondTaskJacoby86(const vec& arg) {
	return -4 * arg[4] * (arg[0] - 2 * arg[5]);
}
double secondTaskJacoby87(const vec& arg) {
	return -arg[1] * exp(arg[1] * arg[6] + arg[9]);
}
double secondTaskJacoby88(const vec& arg) {
	return 0;
}
double secondTaskJacoby89(const vec& arg) {
	return 2 * cos(-arg[8] + arg[2]);
}
double secondTaskJacoby810(const vec& arg) {
	return -exp(arg[1] * arg[6] + arg[9]);
}
double secondTaskJacoby91(const vec& arg) {
	return -3;
}
double secondTaskJacoby92(const vec& arg) {
	return -2 * arg[7] * arg[9] * arg[6];
}
double secondTaskJacoby93(const vec& arg) {
	return 0;
}
double secondTaskJacoby94(const vec& arg) {
	return exp(arg[4] + arg[3]);
}
double secondTaskJacoby95(const vec& arg) {
	return exp(arg[4] + arg[3]);
}
double secondTaskJacoby96(const vec& arg) {
	return -7 / (arg[5] * arg[5]);
}
double secondTaskJacoby97(const vec& arg) {
	return -2 * arg[1] * arg[7] * arg[9];
}
double secondTaskJacoby98(const vec& arg) {
	return -2 * arg[1] * arg[6] * arg[9];
}
double secondTaskJacoby99(const vec& arg) {
	return 3;
}
double secondTaskJacoby910(const vec& arg) {
	return -2 * arg[1] * arg[6] * arg[7];
}
double secondTaskJacoby101(const vec& arg) {
	return arg[9];
}
double secondTaskJacoby102(const vec& arg) {
	return arg[8];
}
double secondTaskJacoby103(const vec& arg) {
	return -arg[7];
}
double secondTaskJacoby104(const vec& arg) {
	return cos(arg[3] + arg[4] + arg[5])*arg[6];
}
double secondTaskJacoby105(const vec& arg) {
	return cos(arg[3] + arg[4] + arg[5])*arg[6];
}
double secondTaskJacoby106(const vec& arg) {
	return cos(arg[3] + arg[4] + arg[5])*arg[6];
}
double secondTaskJacoby107(const vec& arg) {
	return sin(arg[3] + arg[4] + arg[5]);
}
double secondTaskJacoby108(const vec& arg) {
	return -arg[2];
}
double secondTaskJacoby109(const vec& arg) {
	return arg[1];
}
double secondTaskJacoby1010(const vec& arg) {
	return arg[0];
}
*/

double secondTask1(const vec& arg) { return (cos(arg[0] * arg[1]) - exp(-3 * arg[2]) + arg[3] * pow(arg[4], 2) - arg[5] - sinh(2 * arg[7]) * arg[8] + 2 * arg[9] + 2.0004339741653854440); }
double secondTask2(const vec& arg) { return (sin(arg[0] * arg[1]) + arg[2] * arg[8] * arg[6] - exp(-arg[9] + arg[5]) + 3 * pow(arg[4], 2) - arg[5] * (arg[7] + 1) + 10.886272036407019994); }
double secondTask3(const vec& arg) { return (arg[0] - arg[1] + arg[2] - arg[3] + arg[4] - arg[5] + arg[6] - arg[7] + arg[8] - arg[9] - 3.1361904761904761904); }
double secondTask4(const vec& arg) { return (2 * cos(-arg[8] + arg[3]) + arg[4] / (arg[2] + arg[0]) - sin(pow(arg[1], 2)) + pow(cos(arg[6] * arg[9]), 2) - arg[7] - 0.1707472705022304757); }
double secondTask5(const vec& arg) { return (sin(arg[4]) + 2 * arg[7] * (arg[2] + arg[0]) - exp(-arg[6] * (-arg[9] + arg[5])) + 2 * cos(arg[1]) - 1 / (arg[3] - arg[8]) - 0.3685896273101277862); }
double secondTask6(const vec& arg) { return (exp(arg[0] - arg[3] - arg[8]) + pow(arg[4], 2) / arg[7] + cos(3 * arg[9] * arg[1]) / 2 - arg[5] * arg[2] + 2.0491086016771875115); }
double secondTask7(const vec& arg) { return (pow(arg[1], 3) * arg[6] - sin(arg[9] / arg[4] + arg[7]) + (arg[0] - arg[5]) * cos(arg[3]) + arg[2] - 0.7380430076202798014); }
double secondTask8(const vec& arg) { return (arg[4] * pow(arg[0] - 2 * arg[5], 2) - 2 * sin(-arg[8] + arg[2]) + 1.5 * arg[3] - exp(arg[1] * arg[6] + arg[9]) + 3.5668321989693809040); }
double secondTask9(const vec& arg) { return (7 / arg[5] + exp(arg[4] + arg[3]) - 2 * arg[1] * arg[7] * arg[9] * arg[6] + 3 * arg[8] - 3 * arg[0] - 8.4394734508383257499); }
double secondTask10(const vec& arg) { return (arg[9] * arg[0] + arg[8] * arg[1] - arg[7] * arg[2] + sin(arg[3] + arg[4] + arg[5]) * arg[6] - 0.78238095238095238096); }
double secondTaskJacoby11(const vec& arg) { return -sin(arg[0] * arg[1]) * arg[1]; }
double secondTaskJacoby12(const vec& arg) { return -sin(arg[0] * arg[1]) * arg[0]; }
double secondTaskJacoby13(const vec& arg) { return 0.3e1 * exp(-(double)(3 * arg[2])); }
double secondTaskJacoby14(const vec& arg) { return arg[4] * arg[4]; }
double secondTaskJacoby15(const vec& arg) { return 2 * arg[3] * arg[4]; }
double secondTaskJacoby16(const vec& arg) { return -1; }
double secondTaskJacoby17(const vec& arg) { return 0; }
double secondTaskJacoby18(const vec& arg) { return -0.2e1 * cosh((double)(2 * arg[7])) * arg[8]; }
double secondTaskJacoby19(const vec& arg) { return -sinh((double)(2 * arg[7])); }
double secondTaskJacoby110(const vec& arg) { return 2; }
double secondTaskJacoby21(const vec& arg) { return cos(arg[0] * arg[1]) * arg[1]; }
double secondTaskJacoby22(const vec& arg) { return cos(arg[0] * arg[1]) * arg[0]; }
double secondTaskJacoby23(const vec& arg) { return arg[8] * arg[6]; }
double secondTaskJacoby24(const vec& arg) { return 0; }
double secondTaskJacoby25(const vec& arg) { return 6 * arg[4]; }
double secondTaskJacoby26(const vec& arg) { return -exp(-arg[9] + arg[5]) - (double)arg[7] - 0.1e1; }
double secondTaskJacoby27(const vec& arg) { return (double)arg[2] * arg[8]; }
double secondTaskJacoby28(const vec& arg) { return -arg[5]; }
double secondTaskJacoby29(const vec& arg) {
	return (double)arg[2] * arg[6];
}
double secondTaskJacoby210(const vec& arg) {
	return exp(-arg[9] + arg[5]);
}
double secondTaskJacoby31(const vec& arg) { return 1; }
double secondTaskJacoby32(const vec& arg) { return -1; }
double secondTaskJacoby33(const vec& arg) { return 1; }
double secondTaskJacoby34(const vec& arg) { return -1; }
double secondTaskJacoby35(const vec& arg) { return 1; }
double secondTaskJacoby36(const vec& arg) { return -1; }
double secondTaskJacoby37(const vec& arg) { return 1; }
double secondTaskJacoby38(const vec& arg) { return -1; }
double secondTaskJacoby39(const vec& arg) { return 1; }
double secondTaskJacoby310(const vec& arg) { return -1; }
double secondTaskJacoby41(const vec& arg) { return -(double)arg[4] * pow((double)arg[2] + arg[0], -0.2e1); }
double secondTaskJacoby42(const vec& arg) { return -0.2e1 * cos(arg[1] * arg[1]) * arg[1]; }
double secondTaskJacoby43(const vec& arg) { return -(double)arg[4] * pow((double)arg[2] + arg[0], -0.2e1); }
double secondTaskJacoby44(const vec& arg) { return -0.2e1 * sin(-arg[8] + (double)arg[3]); }
double secondTaskJacoby45(const vec& arg) { return 0.1e1 / ((double)arg[2] + arg[0]); }
double secondTaskJacoby46(const vec& arg) { return 0; }
double secondTaskJacoby47(const vec& arg) { return -0.2e1 * cos(arg[6] * arg[9]) * sin(arg[6] * arg[9]) * arg[9]; }
double secondTaskJacoby48(const vec& arg) { return -1; }
double secondTaskJacoby49(const vec& arg) { return 0.2e1 * sin(-arg[8] + (double)arg[3]); }
double secondTaskJacoby410(const vec& arg) { return -0.2e1 * cos(arg[6] * arg[9]) * sin(arg[6] * arg[9]) * arg[6]; }
double secondTaskJacoby51(const vec& arg) { return 2 * arg[7]; }
double secondTaskJacoby52(const vec& arg) { return -0.2e1 * sin(arg[1]); }
double secondTaskJacoby53(const vec& arg) { return 2 * arg[7]; }
double secondTaskJacoby54(const vec& arg) { return pow(-arg[8] + (double)arg[3], -0.2e1); }
double secondTaskJacoby55(const vec& arg) { return cos((double)arg[4]); }
double secondTaskJacoby56(const vec& arg) { return arg[6] * exp(-arg[6] * (-arg[9] + arg[5])); }
double secondTaskJacoby57(const vec& arg) { return -(arg[9] - arg[5]) * exp(-arg[6] * (-arg[9] + arg[5])); }
double secondTaskJacoby58(const vec& arg) { return (double)(2 * arg[2]) + 0.2e1 * arg[0]; }
double secondTaskJacoby59(const vec& arg) { return -pow(-arg[8] + (double)arg[3], -0.2e1); }
double secondTaskJacoby510(const vec& arg) { return -arg[6] * exp(-arg[6] * (-arg[9] + arg[5])); }
double secondTaskJacoby61(const vec& arg) { return exp(arg[0] - (double)arg[3] - arg[8]); }
double secondTaskJacoby62(const vec& arg) { return -0.3e1 / 0.2e1 * sin(0.3e1 * arg[9] * arg[1]) * arg[9]; }
double secondTaskJacoby63(const vec& arg) { return -arg[5]; }
double secondTaskJacoby64(const vec& arg) { return -exp(arg[0] - (double)arg[3] - arg[8]); }
double secondTaskJacoby65(const vec& arg) { return 2 * arg[4] / arg[7]; }
double secondTaskJacoby66(const vec& arg) { return -arg[2]; }
double secondTaskJacoby67(const vec& arg) { return 0; }
double secondTaskJacoby68(const vec& arg) { return -arg[4] * arg[4] * pow((double)arg[7], (double)(-2)); }
double secondTaskJacoby69(const vec& arg) { return -exp(arg[0] - (double)arg[3] - arg[8]); }
double secondTaskJacoby610(const vec& arg) { return -0.3e1 / 0.2e1 * sin(0.3e1 * arg[9] * arg[1]) * arg[1]; }
double secondTaskJacoby71(const vec& arg) { return cos((double)arg[3]); }
double secondTaskJacoby72(const vec& arg) { return 0.3e1 * arg[1] * arg[1] * arg[6]; }
double secondTaskJacoby73(const vec& arg) { return 1; }
double secondTaskJacoby74(const vec& arg) { return -(arg[0] - arg[5]) * sin((double)arg[3]); }
double secondTaskJacoby75(const vec& arg) { return cos(arg[9] / (double)arg[4] + (double)arg[7]) * arg[9] * (double)pow((double)arg[4], (double)(-2)); }
double secondTaskJacoby76(const vec& arg) { return -cos((double)arg[3]); }
double secondTaskJacoby77(const vec& arg) { return pow(arg[1], 0.3e1); }
double secondTaskJacoby78(const vec& arg) { return -cos(arg[9] / (double)arg[4] + (double)arg[7]); }
double secondTaskJacoby79(const vec& arg) { return 0; }
double secondTaskJacoby710(const vec& arg) { return -cos(arg[9] / (double)arg[4] + (double)arg[7]) / (double)arg[4]; }
double secondTaskJacoby81(const vec& arg) { return 0.2e1 * (double)arg[4] * (arg[0] - 0.2e1 * arg[5]); }
double secondTaskJacoby82(const vec& arg) { return -arg[6] * exp(arg[1] * arg[6] + arg[9]); }
double secondTaskJacoby83(const vec& arg) { return -0.2e1 * cos(-arg[8] + (double)arg[2]); }
double secondTaskJacoby84(const vec& arg) { return 0.15e1; }
double secondTaskJacoby85(const vec& arg) { return pow(arg[0] - 0.2e1 * arg[5], 0.2e1); }
double secondTaskJacoby86(const vec& arg) { return -0.4e1 * (double)arg[4] * (arg[0] - 0.2e1 * arg[5]); }
double secondTaskJacoby87(const vec& arg) { return -arg[1] * exp(arg[1] * arg[6] + arg[9]); }
double secondTaskJacoby88(const vec& arg) { return 0; }
double secondTaskJacoby89(const vec& arg) { return 0.2e1 * cos(-arg[8] + (double)arg[2]); }
double secondTaskJacoby810(const vec& arg) { return -exp(arg[1] * arg[6] + arg[9]); }
double secondTaskJacoby91(const vec& arg) { return -3; }
double secondTaskJacoby92(const vec& arg) { return -0.2e1 * (double)arg[7] * arg[9] * arg[6]; }
double secondTaskJacoby93(const vec& arg) { return 0; }
double secondTaskJacoby94(const vec& arg) { return exp((double)(arg[4] + arg[3])); }
double secondTaskJacoby95(const vec& arg) { return exp((double)(arg[4] + arg[3])); }
double secondTaskJacoby96(const vec& arg) { return -0.7e1 * pow(arg[5], -0.2e1); }
double secondTaskJacoby97(const vec& arg) { return -0.2e1 * arg[1] * (double)arg[7] * arg[9]; }
double secondTaskJacoby98(const vec& arg) { return -0.2e1 * arg[1] * arg[9] * arg[6]; }
double secondTaskJacoby99(const vec& arg) { return 3; }
double secondTaskJacoby910(const vec& arg) { return -0.2e1 * arg[1] * (double)arg[7] * arg[6]; }
double secondTaskJacoby101(const vec& arg) { return arg[9]; }
double secondTaskJacoby102(const vec& arg) { return arg[8]; }
double secondTaskJacoby103(const vec& arg) { return -arg[7]; }
double secondTaskJacoby104(const vec& arg) { return cos((double)arg[3] + (double)arg[4] + arg[5]) * arg[6]; }
double secondTaskJacoby105(const vec& arg) { return cos((double)arg[3] + (double)arg[4] + arg[5]) * arg[6]; }
double secondTaskJacoby106(const vec& arg) { return cos((double)arg[3] + (double)arg[4] + arg[5]) * arg[6]; }
double secondTaskJacoby107(const vec& arg) { return sin((double)arg[3] + (double)arg[4] + arg[5]); }
double secondTaskJacoby108(const vec& arg) { return -arg[2]; }
double secondTaskJacoby109(const vec& arg) { return arg[1]; }
double secondTaskJacoby1010(const vec& arg) { return arg[0]; }




vector<funn> secondTaskVec = { secondTask1, secondTask2, secondTask3, secondTask4, secondTask5, secondTask6, secondTask7, secondTask8, secondTask9, secondTask10 };
vector< vector<funn> > secondTaskJacoby = { { secondTaskJacoby11, secondTaskJacoby12, secondTaskJacoby13, secondTaskJacoby14, secondTaskJacoby15, secondTaskJacoby16, secondTaskJacoby17, secondTaskJacoby18, secondTaskJacoby19, secondTaskJacoby110, },
											{ secondTaskJacoby21, secondTaskJacoby22, secondTaskJacoby23, secondTaskJacoby24, secondTaskJacoby25, secondTaskJacoby26, secondTaskJacoby27, secondTaskJacoby28, secondTaskJacoby29, secondTaskJacoby210, },
											{ secondTaskJacoby31, secondTaskJacoby32, secondTaskJacoby33, secondTaskJacoby34, secondTaskJacoby35, secondTaskJacoby36, secondTaskJacoby37, secondTaskJacoby38, secondTaskJacoby39, secondTaskJacoby310, },
											{ secondTaskJacoby41, secondTaskJacoby42, secondTaskJacoby43, secondTaskJacoby44, secondTaskJacoby45, secondTaskJacoby46, secondTaskJacoby47, secondTaskJacoby48, secondTaskJacoby49, secondTaskJacoby410, },
											{ secondTaskJacoby51, secondTaskJacoby52, secondTaskJacoby53, secondTaskJacoby54, secondTaskJacoby55, secondTaskJacoby56, secondTaskJacoby57, secondTaskJacoby58, secondTaskJacoby59, secondTaskJacoby510, },
											{ secondTaskJacoby61, secondTaskJacoby62, secondTaskJacoby63, secondTaskJacoby64, secondTaskJacoby65, secondTaskJacoby66, secondTaskJacoby67, secondTaskJacoby68, secondTaskJacoby69, secondTaskJacoby610, },
											{ secondTaskJacoby71, secondTaskJacoby72, secondTaskJacoby73, secondTaskJacoby74, secondTaskJacoby75, secondTaskJacoby76, secondTaskJacoby77, secondTaskJacoby78, secondTaskJacoby79, secondTaskJacoby710, },
											{ secondTaskJacoby81, secondTaskJacoby82, secondTaskJacoby83, secondTaskJacoby84, secondTaskJacoby85, secondTaskJacoby86, secondTaskJacoby87, secondTaskJacoby88, secondTaskJacoby89, secondTaskJacoby810, },
											{ secondTaskJacoby91, secondTaskJacoby92, secondTaskJacoby93, secondTaskJacoby94, secondTaskJacoby95, secondTaskJacoby96, secondTaskJacoby97, secondTaskJacoby98, secondTaskJacoby99, secondTaskJacoby910, },
											{ secondTaskJacoby101, secondTaskJacoby102, secondTaskJacoby103, secondTaskJacoby104, secondTaskJacoby105, secondTaskJacoby106, secondTaskJacoby107, secondTaskJacoby108, secondTaskJacoby109, secondTaskJacoby1010, }

										   };

#endif