// layer.cpp            V.Rao, H.Rao
// compile for doubleing point hardware if available

#include <cassert>

#include "layer.h"


inline double squash(double input)
// squashing function
// use sigmoid -- can customize to something
// else if desired; can add a bias term too
//
{
    if(input < -50)
    {
        return 0.0;
    }
    else if(input > 50)
    {
        return 1.0;
    }
    else
    {
        return (double)(1 / (1 + exp(-(double)input)));
    }
}

// the next function is needed for Turbo C++
// and Borland C++ to link in the appropriate
// functions for fscanf doubleing point formats:

/*
static void force_fpf()
{
	double x, *y;
	y=&x;
	x=*y;
}
*/

// -----------------------------------------
//  input layer                              
//------------------------------------------
input_layer::input_layer(int i, int o)
{
    num_inputs = i;
    num_outputs = o;
    outputs = NULL;
    outputs = new double[num_outputs];

    assert(outputs);
    orig_outputs = NULL;
    orig_outputs = new double[num_outputs];

    assert(orig_outputs);
    if((outputs == 0) || (orig_outputs == 0))
    {
        cout << "not enough memory\n";
        cout << "choose a smaller architecture\n";
        exit(1);
    }

}

input_layer::~input_layer()
{
    delete[]outputs;
    delete[]orig_outputs;
}


void input_layer::calc_out()
{

    int i;

    for(i = 0; i < num_outputs; i++)
        outputs[i] = orig_outputs[i];

}

// -----------------------------------------
//                              output layer                              
//------------------------------------------ 


output_layer::output_layer(int ins, int outs)
{
    num_inputs = ins;
    num_outputs = outs;
    weights = NULL;
    weights = new double[num_inputs * num_outputs];

    assert(weights);
    outputs = NULL;
    outputs = new double[num_outputs];

    assert(outputs);
}

output_layer::~output_layer()
{
// some compilers may require the array
// size in the delete statement; those
// conforming to Ansi C++ will not
    delete[]weights;
    delete[]outputs;

}


void output_layer::calc_out()
{

    int i, j, k;
    double accumulator = 0.0;


    for(j = 0; j < num_outputs; j++)
    {

        for(i = 0; i < num_inputs; i++)

        {
            k = i * num_outputs;
            if(weights[k + j] * weights[k + j] > 1000000.0)
            {
                cout << "weights are blowing up\n";
                cout << "try a smaller learning constant\n";
                cout << "e.g. beta=0.02    aborting...\n";
                exit(1);
            }
            outputs[j] = weights[k + j] * (*(inputs + i));
            accumulator += outputs[j];
        }
        // use the sigmoid squash function
        outputs[j] = squash(accumulator);
        accumulator = 0;
    }

}

void output_layer::read_weights(int layer_no, FILE * weights_file_ptr)
{
    int i, j, k;


// assume file is already open and ready for
// reading

// look for the prepended layer_no 
// format:
//              layer_no        weight[0,0] weight[0,1] ...
//              layer_no        weight[1,0] weight[1,1] ...
//              ...
    while(1)

    {

        fscanf(weights_file_ptr, "%i", &j);
        if((j == layer_no) || (feof(weights_file_ptr)))
            break;
        else
        {
            while(fgetc(weights_file_ptr) != '\n')
            {;
            }                   // get rest of line
        }
    }

    if(!(feof(weights_file_ptr)))
    {
        // continue getting first line
        i = 0;
        for(j = 0; j < num_outputs; j++)
        {

            fscanf(weights_file_ptr, "%lf", &weights[j]);       // i*num_outputs = 0
        }
        fscanf(weights_file_ptr, "\n");




        // now get the other lines
        for(i = 1; i < num_inputs; i++)
        {
            fscanf(weights_file_ptr, "%i", &layer_no);
            k = i * num_outputs;
            for(j = 0; j < num_outputs; j++)
            {
                fscanf(weights_file_ptr, "%lf", &weights[k + j]);
            }

        }
        fscanf(weights_file_ptr, "\n");
    }


    else
        cout << "end of file reached\n";

}


// -----------------------------------------
//                              middle layer                              
//------------------------------------------


middle_layer::middle_layer(int i, int o):output_layer(i, o)
{

}

middle_layer::~middle_layer()
{
    delete[]weights;
    delete[]outputs;
}


network::network()
{
    // DW 8/24/04: position is never used, so I commented it out
    //position = 0L;
}

network::~network()
{
    delete[]buffer;
}

void network::set_training(const unsigned &value)
{
    training = value;
}

unsigned network::get_training_value()
{
    return training;
}


void network::set_up_network()
{
    int i, j, k;

//-------------------------------------------------------       
// Construct the layers                                                 
//
//-------------------------------------------------------        


    layer_ptr[0] = NULL;
    layer_ptr[0] = new input_layer(0, layer_size[0]);
    assert(layer_ptr[0]);

    for(i = 0; i < (number_of_layers - 1); i++)
    {
        layer_ptr[i + 1] = NULL;
        layer_ptr[i + 1] =
            new middle_layer(layer_size[i], layer_size[i + 1]);
        assert(layer_ptr[i + 1]);
    }
    layer_ptr[number_of_layers - 1] = NULL;
    layer_ptr[number_of_layers - 1] = new
        output_layer(layer_size[number_of_layers - 2],
        layer_size[number_of_layers - 1]);
    assert(layer_ptr[number_of_layers - 1]);
    for(i = 0; i < (number_of_layers - 1); i++)
    {
        if(layer_ptr[i] == 0)
        {
            cout << "insufficient memory\n";
            cout << "use a smaller architecture\n";
            exit(1);
        }
    }

//-------------------------------------------------------       
// Connect the layers
//
//-------------------------------------------------------        
// set inputs to previous layer outputs for all layers,
//      except the input layer

    for(i = 1; i < number_of_layers; i++)
        layer_ptr[i]->inputs = layer_ptr[i - 1]->outputs;



// define the IObuffer that caches data from
// the datafile
    i = layer_ptr[0]->num_outputs;      // inputs
    j = layer_ptr[number_of_layers - 1]->num_outputs;   //outputs
    k = MAX_VECTORS;
    buffer = NULL;
    buffer = new double[(i + j) * k];

    if(buffer == 0)
    {
        cout << "insufficient memory for buffer\n";
        exit(1);
    }
}

void network::read_weights(FILE * weights_file_ptr)
{
    int i;

    for(i = 1; i < number_of_layers; i++)
        ((output_layer *) layer_ptr[i])->read_weights(i, weights_file_ptr);
}


double *network::get_outputs()
{
    return layer_ptr[number_of_layers - 1]->outputs;
}

void network::set_up_pattern(int buffer_index)
{
// read one vector into the network
    int i, k;
    int ins;

    ins = layer_ptr[0]->num_outputs;
    k = buffer_index * ins;

    for(i = 0; i < ins; i++)
        ((input_layer *) layer_ptr[0])->orig_outputs[i] = buffer[k + i];

}

void network::forward_prop()
{
    int i;

    for(i = 0; i < number_of_layers; i++)
    {
        layer_ptr[i]->calc_out();       //polymorphic
        // function
    }
}

layer::~layer()
{
}
