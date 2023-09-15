#pragma once
#ifndef _COMMON_2D_H_
#define _COMMON_2D_H_

// ---------
// MJ include files
#include <vector>
#include <map>
#include <string>
#include "structures.h"


using std::string;
// ---------

namespace fem_2D
{

	class calculation
	{
		public:
			static calculation& getInstance()
			{
				static calculation    instance; // Guaranteed to be destroyed.
									            // Instantiated on first use.
				return instance;
			}
			
		private:
			calculation() 
			{				 
			};                   // Constructor? (the {} brackets) are needed here.
			

			// C++ 11
			// =======
			// We can use the better technique of deleting the methods
			// we don't want.
		public:
			calculation(calculation const&)     = delete;
			void operator=(calculation const&)  = delete;

			// Note: Scott Meyers mentions in his Effective Modern
			//       C++ book, that deleted functions should generally
			//       be public as it results in better error messages
			//       due to the compilers behavior to check accessibility
			//       before deleted status

			
			//WHAT IS BELOW CANNOT BE ACCESED WITHOUT GET INSTANCE
			conf_stru c;

// 			double energy(struct cella_stru& );
// 			struct cella_stru force(struct cella_stru&, struct base_stru&);
	
			void setSize(int n,int nx,int ny)
			{
				c.p.resize(n);
				c.pfix.resize(n);
				c.disorder_x.resize(n);

				c.disp.resize(n);
				c.full.resize(n);

				c.energy_local.resize(n);

				c.gr.resize(n);

				c.grp1.resize(0);
				c.grp1.resize(n);
				c.bc_cond.resize(n);
				c.current_metrics.resize(n);
				c.stress.resize(n);


				//if not initiated properly, it can cause a bug for c.disorder = 0

				for (int i = 0; i < n; ++i){
                    c.grp1[i].resize(12);
                    c.bc_cond[i].resize(4);
                    c.energy_local[i].resize(4);
                    c.current_metrics[i].resize(4);
                    c.stress[i].resize(4);

                }

 

			}



			
	};


}




#endif //


