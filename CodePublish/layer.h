// layer.h              V.Rao, H. Rao
// header file for the layer class heirarchy and
// the network class

#ifndef LAYER_H
#define LAYER_H

#define MAX_LAYERS	5
#define MAX_VECTORS	4


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;


class network;
class Kohonen_network;

class layer
{

protected:

    int num_inputs;
    int num_outputs;
    double *outputs;            // pointer to array of outputs
    double *inputs;             // pointer to array of inputs, which
    // are outputs of some other layer

    friend class network;
    friend class Kohonen_network;     // update for Kohonen model

public:

    virtual void calc_out() = 0;

    virtual ~layer();
};



class input_layer:public layer
{

private:

    double *orig_outputs;

public:

    input_layer(int, int);

    // DW 8/24/04: changed destructor to be virtual
    virtual ~ input_layer();

    virtual void calc_out();

    friend class network;
};

class middle_layer;

class output_layer:public layer
{
protected:

    double *weights;

    friend class network;

public:


    output_layer(int, int);

    // DW 8/24/04: changed destructor to be virtual
    virtual ~ output_layer();

    virtual void calc_out();
    void read_weights(int, FILE *);
};




class middle_layer:public output_layer
{

private:

public:
    middle_layer(int, int);
    ~middle_layer();
};


class network
{
public:

    int number_of_layers;
    int layer_size[MAX_LAYERS];
    double *buffer;

private:

    layer * layer_ptr[MAX_LAYERS];
    
    // DW 8/24/04:  position is never used, so I commented it out
    //fpos_t position;
    
    unsigned training;

public:
    network();
    ~network();
    void set_training(const unsigned &);
    unsigned get_training_value();
    void set_up_network();
    void read_weights(FILE *);
    double *get_outputs();
    void forward_prop();
    void set_up_pattern(int);

};

#endif
