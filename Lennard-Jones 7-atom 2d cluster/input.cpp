#include"input.hpp"
//#include"algorithm_natoms.hpp"

template <class T>
void input::read_item(istream& from, const char *_titre, T * val)
{
  char titre[80];
  char c;

  from.get(titre,80,':');
  from.get(c);
  if(strcmp(titre,_titre))
    {
      cerr<<" The item is titled'"<<titre<<"' while we waited '"<<_titre<<"'."<<endl;
      exit(1);
    }
  from>>*val;
  from.get(titre,80,'\n');
  from.get(c);
}

// void input::read_item_string (istream& from,char *_titre, char *val)
// {
//   char titre[80];
//   char c;
//
//   from.get(titre,80,':');
//   from.get(c);
//
//   if(strcmp(titre,_titre))
//     {
//       cerr<<"L'item a pour titre '"<<titre<<"' alors qu'on attendait '"<<_titre<<"'."<<endl;
//       exit(1);
//     }
//
//   from.get(c);
//   from.get(val,80,'\n');
//   from.get(c);
// }




void input::read(int & Nb_para_, int & N_Coarse_steps_, int & N_Fine_steps_, double& t_step_, int& type_algo_, double& beta_, double& gamma_)
{
  ifstream from("input_file");

  read_item<int>(from,    "Number of para iterations   ",&Nb_para_);
  read_item<int>(from,    "Number of coarse intervals  ",&N_Coarse_steps_);
  read_item<int>(from,    "Number of fine steps        ",&N_Fine_steps_);
  read_item<double>(from, "Time step                   ",&t_step_);
  read_item<int>(from,    "Algorithm type              ",&type_algo_);
  read_item<double>(from, "Beta                        ",&beta_);
  read_item<double>(from, "Gamma                       ",&gamma_);
  // read_item<int>(from,    "Number of atoms             ",&natoms);
  // TODO add new parameters to be read from inp file here

}

void input::load(void)
{
   read(Nb_para,N_Coarse_steps,N_Fine_steps,t_step,type_algo,beta,gamma);
}
