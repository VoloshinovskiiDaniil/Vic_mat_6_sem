#include "func.h"
#include<fstream>

double u_0(double x)
{
	double s = 1 / ((1 + (x - 4.5) * (x - 4.5))) - 4 / 85;
		return s;
}

int main()
{
	const double L = 20;
	const double h = 0.5;
	const double T = 18;
	const double a = 0.3;
	const double Co[3]{ 0.6, 1.0, 1.01 };
	const double tau[3]{ Co[0] * h / a, Co[1] * h / a, Co[2] * h / a };

	////////////////////////////////////////////////////////////
	std::ofstream data_L_angle("data_L_angle_1.01.txt");
	data_L_angle.precision(16);
	std::ofstream data_LW("data_LW_1.01.txt");
	data_LW.precision(16);

	mesh_with_time res_L_angle = Prototype(L, h, u_0);
	for (double t = 0; t < T; t += tau[2])
	{
		data_L_angle << res_L_angle.t << "		";
		for (int i = 0; i < res_L_angle.layer.size(); i++)
		{
			data_L_angle << res_L_angle.layer[i].u << "		";
		}
		data_L_angle << std::endl;
		res_L_angle = L_angle_method(res_L_angle, a, tau[2], h);
	}

	mesh_with_time res_LW = Prototype(L, h, u_0);
	for (double t = 0; t < T; t += tau[2])
	{
		data_LW << res_LW.t << "		";
		for (int i = 0; i < res_LW.layer.size(); i++)
		{
			data_LW << res_LW.layer[i].u << "		";
		}
		data_LW << std::endl;
		res_LW = LW_method(res_LW, a, tau[2], h);
	}

	data_L_angle.close();
	data_LW.close();
	///////////////////////////////////////////////////////////////

	return 0;
}