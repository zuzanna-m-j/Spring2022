from imports import*

def length_geom(system, mc, phi, scale = 0.99, tot_part_vol = 1.0):
        # Calculate initial system packing fraction
        V_init = system.box.get_volume()
        Lx_init = system.box.Lx
        Ly_init = system.box.Ly
        Lz_init = system.box.Lz
        phi_init = tot_part_vol/V_init
        if hoomd.comm.get_rank() == 0:
            print("Beginning fast compression routine")
            print("Current packing fraction is {}".format(phi_init))
        # this is equivalent to V_init/V_final
        alpha = pow(phi_init/phi, 1./2.)

        # run for a hot second to enable overlap checking
        hoomd.run(1)

        # Test whether to expand or compress
        if alpha > 1:
            # Expand the box
            if hoomd.comm.get_rank() == 0: print(" Expanding the box to {}!".format(phi))
            hoomd.update.box_resize(system.box.Lx*alpha, system.box.Ly*alpha,
                                    system.box.Lz*alpha, system.box.xy, system.box.xz,
                                    system.box.yz, period=None)
        else:
            # Compress the box
            if hoomd.comm.get_rank() == 0: print(" Shrinking the box!")
            curr_phi = phi_init
            while curr_phi < phi and hoomd.get_step() <= 1e6:
                if hoomd.comm.get_rank() == 0: print("Current density is {}".format(curr_phi))
                hoomd.update.box_resize(system.box.Lx*scale, system.box.Ly*scale,
                                        system.box.Lz*scale, system.box.xy,
                                        system.box.xz, system.box.yz, period = None)
                overlaps = mc.count_overlaps()
                while overlaps > 0:
                    hoomd.run(100, quiet = True)
                    overlaps = mc.count_overlaps()
                    if hoomd.comm.get_rank() == 0: print(overlaps, end=' ')
                    sys.stdout.flush()
                print();
                curr_phi = tot_part_vol/system.box.get_volume()

            # update to exactly phi to make things kosher.
            hoomd.update.box_resize(Lx = Lx_init*alpha,
                                    Ly = Ly_init*alpha,
                                    Lz = Lz_init*alpha,
                                    xy = system.box.xy, xz = system.box.xz, yz = system.box.yz, period = None)
            hoomd.run(1)
            overlaps = mc.count_overlaps()
            assert(overlaps == 0)
            if hoomd.comm.get_rank() == 0: print("Compression finished! Current density calculated to be {}.".format(tot_part_vol/system.box.get_volume()))


file = 'hexagonal-lattice-base'
seed = random.randint(1,1e7)
vAs = [(-0.21799013664766403, -0.5023571793344106),(0.43806861572454686, 0.25235264134888885),(-0.21996612497820883, 0.2506258175372971)]
vBs = [(-0.5845343621823524, -0.12771177916537713),(0.4002733778408708, -0.30136003049648535),(0.3889689778694574, 0.26161490951675265),(-0.2689354095542138, 0.24840434931275862)]
vCs = [(0.37900872267547725, -0.35788779037693047),(0.05343979731237958, 0.587630519777832),(-0.28050392925997647, 0.20347157435381596),(-0.2790260172522079, -0.3596149095167526)]
vDs = [(-0.5592972649748391, 0.48861096202561094),(-0.22549069130706167, -0.45403063019074075),(0.4325039292599765, -0.446561949111819),(0.4239575873822218, 0.30637513866397204)]


vAu = [(-0.5220179221316321, -0.16503966450482788),(0.47798207786836777, -0.16503966450482788),(0.04496937597614856, 0.3304527187789308)]
vBu = [(0.4740585927875667, -0.36545839924993373),(-0.025942184562527393, 0.5005665557307537),(-0.44993932042913276, 0.13003315464721044),(-0.01692617377975303, -0.3654588399612012)]
vCu = [(0.021942184562527278, 0.5210690568864235),(-0.4780585927875667, -0.34495589809426397),(0.03095664064471082, -0.34495635498980115),(0.4549544416972373, 0.025576284927716575)]
vDu = [(0.014417099295431987, 0.7420396645048278),(-0.4931218469408215, -0.11958915236087647),(-0.06444990268382092, -0.6188416597147626),(0.506840070215735, -0.1283163647786042)]

vAd = [(0.5220179189906603, 0.16503999736020913),(-0.4779820810084331, 0.16503865095000192),(-0.04496871198060394, -0.3304531493205859)]
vBd = [(-0.4740585897588976, 0.36545668411237353),(0.02594335361558253, -0.5005675976613787),(0.4499399905918503, -0.1300336257040997),(0.01692617680738384, 0.36545778589054234)]
vCd = [(-0.014941021565904156, -0.5570686822378279),(0.4850585897588976, 0.30895694594822487),(-0.023956643673533673, 0.30895671750045606),(-0.4479539458367474, -0.06157649329169301)]
vDd = [(-0.023414729954361082, -0.7730399973602091),(0.4841230561755986, 0.08858950286033201),(0.055450439720314625, 0.5878414330454843),(-0.5158388727304595, 0.09731536891911975)]

# defining area of vertices
area1 = hoomd.dem.utils.area(vertices=vAu, factor=1.0)
area2 = hoomd.dem.utils.area(vertices=vBu, factor=1.0)
area3 = hoomd.dem.utils.area(vertices=vCu, factor=1.0)
area4 = hoomd.dem.utils.area(vertices=vDu, factor=1.0)

######
hoomd.context.initialize("--mode=cpu")

#to establish phi_init, a square unit cell should have the following lattice constant:
#a = pow(1./phi_init, 1./2.)

#system = hoomd.init.create_lattice(unitcell=hoomd.lattice.sq(a=a), n=10)

# c = 0.72
# d = 0.7
# positions =[]
# positions.append((c - 0.451 , d + 0.419, 0))
# positions.append((c + 0.269 , d + 0.419, 0))
# positions.append((c - 0.391 , d - 0.321, 0))
# positions.append((c + 0.239 , d - 0.221, 0))
# uc = hoomd.lattice.unitcell(N = 4,
#                             a1 = [1.4,0,0],
#                             a2 = [0,1.4,0],
#                             a3 = [0,0,1.0],
#                             dimensions = 2,
#                             position = positions,
#                             type_name = ['As', 'Bs', 'Cs', 'Ds']);
#
# Set the log file for the data
# I initialize parameters for the Frenkel-Ladd calculation




sx = 1.9070982572923945/2 + 0.15
ux = 1.9070982572923945/2 + 0.15
sy = 1.7620793290096557 + 0.075

# x: 1.9070982572923945
# y: 1.7620793290096557

positions =[]
positions.append((0.027, -0.716, 0))
positions.append((0.552, -0.511, 0))
positions.append((-0.548, -0.536,0))
positions.append((-0.015, 0.139, 0))
positions.append((-0.027 + sx , 0.716, 0))
positions.append((-0.552 + sx, 0.511, 0))
positions.append((0.541 + sx , 0.572, 0))
positions.append((0.024 + sx, -0.108, 0))

positions.append((0.027 + ux, -0.716 + sy, 0))
positions.append((0.552 + ux, -0.511 + sy, 0))
positions.append((-0.548 + ux, -0.536 + sy,0))
positions.append((-0.015 + ux, 0.139 + sy, 0))
positions.append((-0.027 + sx + ux , 0.716 + sy, 0))
positions.append((-0.552 + sx + ux, 0.511 + sy, 0))
positions.append((0.541 + sx + ux, 0.572 + sy, 0))
positions.append((0.024 + sx + ux, -0.108 + sy, 0))

uc = hoomd.lattice.unitcell(N = 16,
                            a1 = [2.2,0,0],
                            a2 = [0,3.67,0],
                            a3 = [0,0,1.0],
                            dimensions = 2,
                            position = positions,
                            type_name = ['Au', 'Bu', 'Cu', 'Du', 'Ad', 'Bd', 'Cd', 'Dd','Au', 'Bu', 'Cu', 'Du', 'Ad', 'Bd', 'Cd', 'Dd']);


system = hoomd.init.create_lattice(unitcell=uc, n=[14,7])

#defining monte carlo simulation

mc = hoomd.hpmc.integrate.convex_polygon(d=0.2, a=0.2, seed=seed)


# mc.shape_param.set('As', vertices=vAs)
# mc.shape_param.set('Bs', vertices=vBs)
# mc.shape_param.set('Cs', vertices=vCs)
# mc.shape_param.set('Ds', vertices=vDs)

mc.shape_param.set('Au', vertices=vAu)
mc.shape_param.set('Bu', vertices=vBu)
mc.shape_param.set('Cu', vertices=vCu)
mc.shape_param.set('Du', vertices=vDu)
mc.shape_param.set('Ad', vertices=vAd)
mc.shape_param.set('Bd', vertices=vBd)
mc.shape_param.set('Cd', vertices=vCd)
mc.shape_param.set('Dd', vertices=vDd)

# gsd = hoomd.dump.gsd(file + ".gsd",
#                    period=100,
#                    group=hoomd.group.all(),
#                    overwrite=True,)

#gsd.dump_shape(mc)
area1 = hoomd.dem.utils.area(vertices=vAs, factor=1.0)
area2 = hoomd.dem.utils.area(vertices=vBs, factor=1.0)
area3 = hoomd.dem.utils.area(vertices=vCs, factor=1.0)
area4 = hoomd.dem.utils.area(vertices=vDs, factor=1.0)
vol  = (area4 + area3 + area2 + area1)*14*7*4
# hoomd.analyze.log(filename=file + '.txt', quantities=['volume', 'N'], period=10000, overwrite=True)
# hoomd.hpmc.analyze.sdf(mc=mc, filename=file + '.dat', xmax=0.02, dx=1e-4, navg=100, period=100, overwrite=True)
# hoomd.run(1)
# length_geom(system, mc, 95, scale = 0.9995, tot_part_vol = vol)
#
# hoomd.run(100)

# r0 = [ ]
# for i in range(len(system.particles)):
#     p = system.particles.pdata.getPosition(i)
#     r0.append([p.x,p.y,p.z])
# q0 = [(1, 0, 0, 0)] * len(system.particles)
#
# fl = hoomd.hpmc.field.frenkel_ladd_energy(mc=mc, ln_gamma=0.0, q_factor=1.0, r0=r0, q0=q0, drift_period=1)
# log_fl = hoomd.analyze.log(filename= file + '.fl',quantities=["lattice_energy","lattice_translational_spring_constant","lattice_energy_pp_avg","lattice_energy_pp_sigma","lattice_num_samples"],period=1)
# lambda_max = 3000
# ks = np.logspace(np.log10(0.00001),np.log10(lambda_max), 200)
# for k in ks:
#     fl.set_params(ln_gamma=np.log(k), q_factor=1.0)
#     fl.reset_statistics()
#     hoomd.run(1000)


import matplotlib.pyplot as plt
import numpy as np
import math

def extrapolate(s, dx, xmax, degree=5):
  s = s[:int(math.ceil(xmax/dx))]
  x = np.arange(0,xmax,dx)
  x += dx/2
  x = x[:len(s)]
  p = np.polyfit(x, s, degree)
  return np.polyval(p, 0.0)

#full_volume_data = np.loadtxt(file + '.vol',skiprows=1)
sdf_data = np.loadtxt(file + '.dat')
sdf_t = sdf_data[:,0]
sdf_data = sdf_data[:,1:]
volume_data = []
t_steps = list(map(int,sdf_t))

# for vol in full_volume_data:
#     if int(vol[0]) in t_steps:
#         volume_data.append(float(vol[1]))
#
# N = full_volume_data[0,2]
volume_data = system.box.get_volume()
n = len(sdf_data)
dim = 2
xmax = 0.02
dx = 1e-4


s0 = []
for i in range(n):
    s0.append(extrapolate(sdf_data[i], dx, xmax, degree = 5))
s0 = np.array(s0)
N = 16*7*14
n_density = N/volume_data
betaP = n_density*(1.0 + s0/(2*dim))

with open("hex-expand.data",'w') as pv:
    pv.writelines('timestep:    volume:     betaP:\n')
    for i in range(len(t_steps)):
        pv.writelines(str(t_steps[i]) + ' ' + str(volume_data) + '  ' + str(betaP[i]) + '\n')

plt.plot(t_steps,betaP,'--')
plt.show()
