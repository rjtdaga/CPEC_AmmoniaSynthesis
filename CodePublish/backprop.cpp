// backprop.cpp		V. Rao, H. Rao
#include "layer.h"                      
#ifndef NEURAL_NET
#define NEURAL_NET

#define WEIGHTS_FILE "weights.dat"
/* Modification EH : Entry into function */
// LayerSize = [Input Nodes, Intermediate Nodes Layer 1, ..., Output Nodes]
	
// create a network object
static network backp;

void InitializeNetwork(int NumberOfLayers, int LayerSize[])
{
	FILE * weights_file_ptr = NULL;
    
	// enter the training mode : 1=training on     0=training off
	backp.set_training((unsigned)0);
	
    // get layer information
	backp.number_of_layers = NumberOfLayers;
	for(int i = 0; i < NumberOfLayers; ++i)
    {
        backp.layer_size[i] = LayerSize[i];
	}
    
    // set up the network connections
	backp.set_up_network();
    
    // initialize the weights
	// read in the weight matrix defined by a
	// prior run of the backpropagation simulator
	// with training on
    weights_file_ptr=fopen(WEIGHTS_FILE,"r");
	if(weights_file_ptr ==NULL)
	{
		//cout << "\nNeural Network NOT Initialized (weights are absent)"
      //      << endl;
		return; 
	}
	backp.read_weights(weights_file_ptr);
	fclose(weights_file_ptr);
}

double* NeuralNetwork(double *buffer)
{
	backp.buffer = buffer;
	backp.set_up_pattern(0);
	backp.forward_prop();
	double* results = backp.get_outputs();
	return results;
}
#endif
