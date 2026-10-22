function [state, m, name_origin_list, time_start, time_end, top_five_pop_flow_nodes, time_step_size , zeta_vec, zeta,  Nv, gamma, beta_1_vec, phi_0] = county_name_list()

%% variable parameters 

state = 48;

m = 254;

name_origin_list = {'Anderson County', 'Andrews County', 'Angelina County', 'Aransas County', 'Archer County', 'Armstrong County', 'Atascosa County', 'Austin County', 'Bailey County', 'Bandera County', 'Bastrop County', 'Baylor County', 'Bee County', 'Bell County', 'Bexar County', 'Blanco County', 'Borden County', 'Bosque County', 'Bowie County', 'Brazoria County', 'Brazos County', 'Brewster County', 'Briscoe County', 'Brooks County', 'Brown County', 'Burleson County', 'Burnet County', 'Caldwell County', 'Calhoun County', 'Callahan County', 'Cameron County', 'Camp County', 'Carson County', 'Cass County', 'Castro County', 'Chambers County', 'Cherokee County', 'Childress County', 'Clay County', 'Cochran County', 'Coke County', 'Coleman County', 'Collin County', 'Collingsworth County', 'Colorado County', 'Comal County', 'Comanche County', 'Concho County', 'Cooke County', 'Coryell County', 'Cottle County', 'Crane County', 'Crockett County', 'Crosby County', 'Culberson County', 'Dallam County', 'Dallas County', 'Dawson County', 'DeWitt County', 'Deaf Smith County', 'Delta County', 'Denton County', 'Dickens County', 'Dimmit County', 'Donley County', 'Duval County', 'Eastland County', 'Ector County', 'Edwards County', 'El Paso County', 'Ellis County', 'Erath County', 'Falls County', 'Fannin County', 'Fayette County', 'Fisher County', 'Floyd County', 'Foard County', 'Fort Bend County', 'Franklin County', 'Freestone County', 'Frio County', 'Gaines County', 'Galveston County', 'Garza County', 'Gillespie County', 'Glasscock County', 'Goliad County', 'Gonzales County', 'Gray County', 'Grayson County', 'Gregg County', 'Grimes County', 'Guadalupe County', 'Hale County', 'Hall County', 'Hamilton County', 'Hansford County', 'Hardeman County', 'Hardin County', 'Harris County', 'Harrison County', 'Hartley County', 'Haskell County', 'Hays County', 'Hemphill County', 'Henderson County', 'Hidalgo County', 'Hill County', 'Hockley County', 'Hood County', 'Hopkins County', 'Houston County', 'Howard County', 'Hudspeth County', 'Hunt County', 'Hutchinson County', 'Irion County', 'Jack County', 'Jackson County', 'Jasper County', 'Jeff Davis County', 'Jefferson County', 'Jim Hogg County', 'Jim Wells County', 'Johnson County', 'Jones County', 'Karnes County', 'Kaufman County', 'Kendall County', 'Kenedy County', 'Kent County', 'Kerr County', 'Kimble County', 'King County', 'Kinney County', 'Kleberg County', 'Knox County', 'La Salle County', 'Lamar County', 'Lamb County', 'Lampasas County', 'Lavaca County', 'Lee County', 'Leon County', 'Liberty County', 'Limestone County', 'Lipscomb County', 'Live Oak County', 'Llano County', 'Loving County', 'Lubbock County', 'Lynn County', 'Madison County', 'Marion County', 'Martin County', 'Mason County', 'Matagorda County', 'Maverick County', 'McCulloch County', 'McLennan County', 'McMullen County', 'Medina County', 'Menard County', 'Midland County', 'Milam County', 'Mills County', 'Mitchell County', 'Montague County', 'Montgomery County', 'Moore County', 'Morris County', 'Motley County', 'Nacogdoches County', 'Navarro County', 'Newton County', 'Nolan County', 'Nueces County', 'Ochiltree County', 'Oldham County', 'Orange County', 'Palo Pinto County', 'Panola County', 'Parker County', 'Parmer County', 'Pecos County', 'Polk County', 'Potter County', 'Presidio County', 'Rains County', 'Randall County', 'Reagan County', 'Real County', 'Red River County', 'Reeves County', 'Refugio County', 'Roberts County', 'Robertson County', 'Rockwall County', 'Runnels County', 'Rusk County', 'Sabine County', 'San Augustine County', 'San Jacinto County', 'San Patricio County', 'San Saba County', 'Schleicher County', 'Scurry County', 'Shackelford County', 'Shelby County', 'Sherman County', 'Smith County', 'Somervell County', 'Starr County', 'Stephens County', 'Sterling County', 'Stonewall County', 'Sutton County', 'Swisher County', 'Tarrant County', 'Taylor County', 'Terrell County', 'Terry County', 'Throckmorton County', 'Titus County', 'Tom Green County', 'Travis County', 'Trinity County', 'Tyler County', 'Upshur County', 'Upton County', 'Uvalde County', 'Val Verde County', 'Van Zandt County', 'Victoria County', 'Walker County', 'Waller County', 'Ward County', 'Washington County', 'Webb County', 'Wharton County', 'Wheeler County', 'Wichita County', 'Wilbarger County', 'Willacy County', 'Williamson County', 'Wilson County', 'Winkler County', 'Wise County', 'Wood County', 'Yoakum County', 'Young County', 'Zapata County', 'Zavala County'}


time_start = 1;

time_end = 93;

top_five_pop_flow_nodes = [101,57,220,15,227];

%% constant parameters which we may change 

time_step_size = 0.10;                                                     % can compute N_time_steps

zeta_vec = [5, 6, 7, 8, 9, 10];                                            % put the diferent zeta values 
             
zeta = zeta_vec(1) ; % Dimensionless multiplicative factor that modifies the infectivity at the heterogeneous node.

Nv = 3;

gamma = 0.1429;                                                            % 1/7 wekly recovered

beta_1_vec = linspace(0,0.1428,1);                                         % simulate the infection rates

phi_0 = 0.004;                                                             % mean flux value in the first 81 days i.e. 0.006795

end