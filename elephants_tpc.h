#ifndef elephants_tpc_h___
#define elephants_tpc_h___

//             if ( count_years==246 && day_in_year==0 )
            if ( count_years==226 && day_in_year==0 )    // reference: 1987
            {
                EL_TPC_2009 = MyTreePop.getpCanopy();
                for ( int EL_cnt=0; EL_cnt<4; EL_cnt++ ) EL_TPC_size_2009[EL_cnt] = 0;
                for ( int EL_cnt=0; EL_cnt<MyTreePop.getPopSize(); EL_cnt++ )
                {
                    double tmp_height = MyTreePop.getTreeHeight(EL_cnt);
                    if      (                   tmp_height <= 1. ) EL_TPC_size_2009[0]++;
                    else if ( tmp_height >1. && tmp_height <= 3. ) EL_TPC_size_2009[1]++;
                    else if ( tmp_height >3. && tmp_height <= 5. ) EL_TPC_size_2009[2]++;
                    else if ( tmp_height >5.                     ) EL_TPC_size_2009[3]++;
                }
                for ( int EL_cnt=0; EL_cnt<4; EL_cnt++ ) EL_TPC_size_2009[EL_cnt] /= MyTreePop.getPopSize();
                for ( int EL_cnt=0; EL_cnt<4; EL_cnt++ ) cout << EL_TPC_size_2009[EL_cnt] << endl;
            }
            if ( count_years >246 && day_in_year==0 )
            {
//                 if ( ( (MyTreePop.getpCanopy()/EL_TPC_2009-1.)*100. )<(-30.) )
                // TPC: decrease of tree cover by more than 30%
                if ( MyTreePop.getpCanopy() < 0.7*EL_TPC_2009 )
                        EL_TPC_dat << "tc30 "
                            << setw(10) << longitude
                            << setw(10) << latitude
                            << setw( 7) << count_years
                            << setw(14) << EL_TPC_2009
                            << setw(14) << MyTreePop.getpCanopy()
                            << setw(14) << 0
                            << setw(14) << 0
                            << setw(14) << i_scen
                            << setw(14) << with_fire
                            << setw(14) << d_seed
                            << setw(14) << d_param_1
                            << setw(14) << d_param_2
                            << setw(14) << d_param_3
                            << setw(14) << d_param_4
                            << setw(14) << d_param_5
                            << setw(14) << d_param_6
                            << setw(14) << d_param_7
                            << setw(14) << d_param_8
                            << endl;
//                 if ( ( (MyTreePop.getpCanopy()/EL_TPC_2009-1.)*100. )<(-80.) )
                // TPC: decrease of tree cover by more than 80%
                if ( MyTreePop.getpCanopy() < 0.2*EL_TPC_2009 )
                    EL_TPC_dat << "tc80 "
                            << setw(10) << longitude
                            << setw(10) << latitude
                            << setw( 7) << count_years
                            << setw(14) << EL_TPC_2009
                            << setw(14) << MyTreePop.getpCanopy()
                            << setw(14) << 0
                            << setw(14) << 0
                            << setw(14) << i_scen
                            << setw(14) << with_fire
                            << setw(14) << d_seed
                            << setw(14) << d_param_1
                            << setw(14) << d_param_2
                            << setw(14) << d_param_3
                            << setw(14) << d_param_4
                            << setw(14) << d_param_5
                            << setw(14) << d_param_6
                            << setw(14) << d_param_7
                            << setw(14) << d_param_8
                            << endl;

                for ( int EL_cnt=0; EL_cnt<4; EL_cnt++ ) EL_TPC_size_str[EL_cnt] = 0;
                for ( int EL_cnt=0; EL_cnt<MyTreePop.getPopSize(); EL_cnt++ )
                {
                    double tmp_height = MyTreePop.getTreeHeight(EL_cnt);
                    if      (                   tmp_height <= 1. ) EL_TPC_size_str[0]++;
                    else if ( tmp_height >1. && tmp_height <= 3. ) EL_TPC_size_str[1]++;
                    else if ( tmp_height >3. && tmp_height <= 5. ) EL_TPC_size_str[2]++;
                    else if ( tmp_height >5.                     ) EL_TPC_size_str[3]++;
                }

                // TPC: abundance of two size classes > 90%
                if ( (double)(EL_TPC_size_str[0]+EL_TPC_size_str[1])/(double)MyTreePop.getPopSize() > 0.9 ||
                     (double)(EL_TPC_size_str[0]+EL_TPC_size_str[2])/(double)MyTreePop.getPopSize() > 0.9 ||
                     (double)(EL_TPC_size_str[0]+EL_TPC_size_str[3])/(double)MyTreePop.getPopSize() > 0.9 ||
                     (double)(EL_TPC_size_str[1]+EL_TPC_size_str[2])/(double)MyTreePop.getPopSize() > 0.9 ||
                     (double)(EL_TPC_size_str[1]+EL_TPC_size_str[3])/(double)MyTreePop.getPopSize() > 0.9 ||
                     (double)(EL_TPC_size_str[2]+EL_TPC_size_str[3])/(double)MyTreePop.getPopSize() > 0.9 )
                {
                    EL_TPC_dat << "hstr "
                            << setw(10) << longitude
                            << setw(10) << latitude
                            << setw( 7) << count_years
                            << setw(14) << EL_TPC_size_str[0]
                            << setw(14) << EL_TPC_size_str[1]
                            << setw(14) << EL_TPC_size_str[2]
                            << setw(14) << EL_TPC_size_str[3]
                            << setw(14) << i_scen
                            << setw(14) << with_fire
                            << setw(14) << d_seed
                            << setw(14) << d_param_1
                            << setw(14) << d_param_2
                            << setw(14) << d_param_3
                            << setw(14) << d_param_4
                            << setw(14) << d_param_5
                            << setw(14) << d_param_6
                            << setw(14) << d_param_7
                            << setw(14) << d_param_8
                            << endl;
                }

                // TPC: 40% reduction of small trees (<1m)
                if ( ((double)EL_TPC_size_str[0]/(double)MyTreePop.getPopSize()) < 0.6*EL_TPC_size_2009[0] )
                {
                    cout << ((double)EL_TPC_size_str[0]/(double)MyTreePop.getPopSize()) << endl;
                    EL_TPC_dat << "smtr "
                            << setw(10) << longitude
                            << setw(10) << latitude
                            << setw( 7) << count_years
                            << setw(14) << EL_TPC_size_str[0]
                            << setw(14) << EL_TPC_size_str[1]
                            << setw(14) << EL_TPC_size_str[2]
                            << setw(14) << EL_TPC_size_str[3]
                            << setw(14) << i_scen
                            << setw(14) << with_fire
                            << setw(14) << d_seed
                            << setw(14) << d_param_1
                            << setw(14) << d_param_2
                            << setw(14) << d_param_3
                            << setw(14) << d_param_4
                            << setw(14) << d_param_5
                            << setw(14) << d_param_6
                            << setw(14) << d_param_7
                            << setw(14) << d_param_8
                            << endl;
                }
                // TPC: 10% reduction of large trees (>5m)
                if ( ((double)EL_TPC_size_str[3]/(double)MyTreePop.getPopSize()) < 0.9*EL_TPC_size_2009[3] )
                {
                    cout << ((double)EL_TPC_size_str[3]/(double)MyTreePop.getPopSize()) << endl;
                    EL_TPC_dat << "latr "
                            << setw(10) << longitude
                            << setw(10) << latitude
                            << setw( 7) << count_years
                            << setw(14) << EL_TPC_size_str[0]
                            << setw(14) << EL_TPC_size_str[1]
                            << setw(14) << EL_TPC_size_str[2]
                            << setw(14) << EL_TPC_size_str[3]
                            << setw(14) << i_scen
                            << setw(14) << with_fire
                            << setw(14) << d_seed
                            << setw(14) << d_param_1
                            << setw(14) << d_param_2
                            << setw(14) << d_param_3
                            << setw(14) << d_param_4
                            << setw(14) << d_param_5
                            << setw(14) << d_param_6
                            << setw(14) << d_param_7
                            << setw(14) << d_param_8
                            << endl;
                }

            }


#endif

