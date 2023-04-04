#include"initial.hpp"

bool initial() {

	switch (geo)
	{
	case 4:
		ellipsoid();
		//return true;
		break;
	default:
		printf("Non-available geometry!\n");
		return false;
		break;
	}

	FILE* energy;
	FILE* grid;
	energy = fopen("energy.out", "w");
	fprintf(energy,"cycle\tEnergy_diff\tEnergy_ldg\tEnergy_l1\tEnergy_l2\tEnergy_l3\tEnergy_l4\tEnergy_chiral\tEnergy_surf\tEnergy_surf\tEnergy_tot\n");
	fclose(energy);

	int signal;

	grid = fopen("grid.bin", "wb");
	
	for(int l = 0; l < total_points; l++){
		//if(boundary[l] || nboundary[l])	signal = 1;
		if(boundary[l])	signal = 1;
		else if(drop[l]) 	signal = 0;			
		else signal = -1;
		fwrite(&signal, sizeof(int), 1, grid);
	}
	fclose(grid);

	return true;
}