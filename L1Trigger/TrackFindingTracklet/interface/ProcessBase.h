//Base class for processing modules
#ifndef PROCESSBASE_H
#define PROCESSBASE_H

using namespace std;

class ProcessBase{

public:

  ProcessBase(string name, unsigned int iSector){
    name_=name;
    iSector_=iSector;
  }

  virtual ~ProcessBase() { } 

  virtual void addOutput(MemoryBase* memory,string output)=0;

  virtual void addInput(MemoryBase* memory,string input)=0;

  string getName() const {return name_;}

  unsigned int nbits(unsigned int power) {

    if (power==2) return 1;
    if (power==4) return 2;
    if (power==8) return 3;
    if (power==16) return 4;
    if (power==32) return 5;

    cout << "nbits: power = "<<power<<endl;
    assert(0);

    return -1;
    
  }


protected:

  string name_;
  unsigned int iSector_;


};

#endif
