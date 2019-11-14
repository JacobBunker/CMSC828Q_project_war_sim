   clock_t start, end;
   double cpu_time_used1,cpu_time_used2;
   int nneuron=5000,nweights=10000;
   neuralnet *nn1=neuralnet_init(1,1,nneuron,nweights),*nn2=neuralnet_init(1,1,nneuron,nweights);
   neuralnet *nn4=neuralnet_init(1,1,nneuron,nweights),*nn3=neuralnet_init(1,1,nneuron,nweights);

   
   neuralnet_replace(nn3,nn1);
   neuralnet_replace(nn4,nn1);
   
   start = clock();
   neuralnet_replace(nn1,nn2);
   end = clock();
   cpu_time_used1 = ((double) (end - start)) / CLOCKS_PER_SEC;
   start = clock();
   nn3=neuralnet_copy(nn4);
   end = clock();
   cpu_time_used2 = ((double) (end - start)) / CLOCKS_PER_SEC;
   printf("replace %f,copy %f\n",cpu_time_used1,cpu_time_used2);
