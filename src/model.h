#ifndef DEF_MODEL
#define DEF_MODEL

#include <string>
#include <set>
#include <vector>

/* "extern" variables, used in model.cpp, initialized in main.cpp */
extern const std::string natSelFile;   
extern const std::string sexSelFile;
extern const std::string phenoFile;
extern const std::string geneFile;
extern std::string chromoF;

class Model
{
    public:
    Model();                    
    Model(std::string, int, int);
    Model(std::string);
    void writeModelFiles() const;              
    void printModel(std::ostream &) const;
    std::string getModelName() const;
    int getTotalHybPop() const;
    ~Model();                     
    
    private:
    std::string model_name; 
    int n_dist;
    int n_sel_genes;
    std::vector<std::string> list_genes_sel;                          //+1 sex gene
    std::vector<std::string> list_pop_name = {"pop1", "pop2", "hybFoo"};
    std::vector<char> list_pop_lab = {'x', 'y', 'z'};
    void writeSexualSelectionFile() const; // called in writeModelFiles() 
    void writeNaturalSelectionFile() const; // called in writeModelFiles() 
    void writePhenotypeFile() const; // called in writeModelFiles()     
    void writeGeneFile() const;   // called in writeModelFiles() 
};


std::ostream &operator<<( std::ostream &, Model const &);

#endif
