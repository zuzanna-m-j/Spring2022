import numpy
from imports import*

# define logging options

seed = random.randint(0,1e7)
file = sys.argv[0][:-3]
dir = sys.argv[1]
os.system(fr'mkdir {dir}')
with open(fr'{dir}/' + file + '.info', 'w') as info:
    info.write(f'seed: {seed}\n')
    info.writelines(f'filename: {file}\n')
    info.writelines(f'Time: {str(datetime.datetime.now())}')
os.system(fr'cp {file}.py {dir}/script.py')
#sys.stdout = open(fr'{dir}/' + file + '.log', 'w')

# initialize simulation context
hoomd.context.initialize("--mode=cpu")

# Define the vertices for the hexagonal lattice
vAu = [(-0.5220179221316321, -0.16503966450482788),(0.47798207786836777, -0.16503966450482788),(0.04496937597614856, 0.3304527187789308)]
vBu = [(0.4740585927875667, -0.36545839924993373),(-0.025942184562527393, 0.5005665557307537),(-0.44993932042913276, 0.13003315464721044),(-0.01692617377975303, -0.3654588399612012)]
vCu = [(0.021942184562527278, 0.5210690568864235),(-0.4780585927875667, -0.34495589809426397),(0.03095664064471082, -0.34495635498980115),(0.4549544416972373, 0.025576284927716575)]
vDu = [(0.014417099295431987, 0.7420396645048278),(-0.4931218469408215, -0.11958915236087647),(-0.06444990268382092, -0.6188416597147626),(0.506840070215735, -0.1283163647786042)]

vAd = [(0.5220179189906603, 0.16503999736020913),(-0.4779820810084331, 0.16503865095000192),(-0.04496871198060394, -0.3304531493205859)]
vBd = [(-0.4740585897588976, 0.36545668411237353),(0.02594335361558253, -0.5005675976613787),(0.4499399905918503, -0.1300336257040997),(0.01692617680738384, 0.36545778589054234)]
vCd = [(-0.014941021565904156, -0.5570686822378279),(0.4850585897588976, 0.30895694594822487),(-0.023956643673533673, 0.30895671750045606),(-0.4479539458367474, -0.06157649329169301)]
vDd = [(-0.023414729954361082, -0.7730399973602091),(0.4841230561755986, 0.08858950286033201),(0.055450439720314625, 0.5878414330454843),(-0.5158388727304595, 0.09731536891911975)]

# Create a unit cell
sx = 1.9070982572923945/2 + 0.15
ux = 1.9070982572923945/2 + 0.15
sy = 1.7620793290096557 + 0.075

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
                            type_name = ['Au', 'Bu', 'Cu', 'Du', 'Ad', 'Bd', 'Cd', 'Dd','Au', 'Bu', 'Cu', 'Du', 'Ad', 'Bd', 'Cd', 'Dd'])
system = hoomd.init.create_lattice(unitcell=uc, n=[14,7])

# defining area of vertices
area1 = hoomd.dem.utils.area(vertices=vAu, factor=1.0)
area2 = hoomd.dem.utils.area(vertices=vBu, factor=1.0)
area3 = hoomd.dem.utils.area(vertices=vCu, factor=1.0)
area4 = hoomd.dem.utils.area(vertices=vDu, factor=1.0)
tot_vol = (area4 + area2 + area3 + area3)/4

# create the integrator
d_moves = {}
a_moves = {}
for t in system.particles.types:
    d_moves[t] = 0.05
    a_moves[t] = 0.05
mc = hoomd.hpmc.integrate.convex_polygon(d = d_moves, a = a_moves, seed=seed, nselect=4)

# Add vertices to the integrator
mc.shape_param.set('Au', vertices=vAu)
mc.shape_param.set('Bu', vertices=vBu)
mc.shape_param.set('Cu', vertices=vCu)
mc.shape_param.set('Du', vertices=vDu)
mc.shape_param.set('Ad', vertices=vAd)
mc.shape_param.set('Bd', vertices=vBd)
mc.shape_param.set('Cd', vertices=vCd)
mc.shape_param.set('Dd', vertices=vDd)

# MOVE TUNER - tune the acceptacne ratio per particle type
# We expect different acceptance ratios for our particles
tuners = {}
for t in system.particles.types:
    tuners[t] = hpmc.util.tune(mc, tunables=['d', 'a'], type = t, target=0.2, gamma=0.5)

# run the simulation for some time to adjust the move size
for _ in range(10):
    hoomd.run(1000)
    for t in system.particles.types:
        tuners[t].update()

# save the trajectory

gsd = hoomd.dump.gsd(fr'{dir}/' +file + ".gsd",
                   period=1000,
                   group=hoomd.group.all(),
                   overwrite=True,)
gsd.dump_shape(mc)

# Initialize log files
# gather initial pressure statistics

vol_init = hoomd.analyze.log(filename=fr'{dir}/' + file + '-init.vol', quantities=['volume', 'N'], period=100, overwrite=True)
sdf_init = hoomd.hpmc.analyze.sdf(mc=mc, filename=fr'{dir}/' + file + '-init.sdf', xmax=0.001, dx=1e-7, navg=100, period=100, overwrite=True)

print(f'Gather pressure statistics: {hoomd.get_step()}\n')
hoomd.run(5e5)
sdf_init.disable()
vol_init.disable()
print(f"Stop pressure statistics: {hoomd.get_step()}\n")


############################################
# Equilibriate to betaP0

betaP0 = 38.0
boxMC = hpmc.update.boxmc(mc, betaP=betaP0, seed=seed)
boxMC.set_betap(betaP0)
boxMC.volume(delta=system.box.get_volume() * 0.1, weight=1.0)

print(f'Equilibriate to to {betaP0}: {hoomd.get_step()}\n')

sdf_eq = hoomd.hpmc.analyze.sdf(mc=mc, filename=fr'{dir}/' + file + '-eq.sdf', xmax=0.001, dx=1e-7, navg=100, period=100, overwrite=True)
vol_eq = hoomd.analyze.log(filename=fr'{dir}/' + file + '-eq.vol', quantities=['volume', 'N'], period=100, overwrite=True)

tuner = hpmc.util.tune_npt(boxMC, tunables=['dV'], target=0.3, gamma=1.0)
for i in range(10):
    hoomd.run(1e4)
    tuner.update()
sdf_eq.disable()
vol_eq.disable()

print(f'End equilibriation to {betaP0}: {hoomd.get_step()}')


#############################################
# NPT simulation

print(f"NPT simulation: {hoomd.get_step()}")
betaP_logger = hoomd.hpmc.analyze.log(filename='betaP.data', quantities=['volume', 'P', 'N'], period=100)
sdf_npt = hoomd.hpmc.analyze.sdf(mc=mc, filename=fr'{dir}/' + file + '-npt.sdf', xmax=0.005, dx=1e-6, navg=100, period=100, overwrite=True)
for p in numpy.linspace(betaP0,0.0001,100):
    betaP_logger.register_callback('P', lambda timestep: p)
    sdf_eq.disable()
    boxMC.set_betap(p)
    hoomd.run(1e4)
    tuner.update()
    for t in system.particles.types:
        tuners[t].update()
    sdf_npt.enable()
    hoomd.run(1e5)


def extrapolate(s, dx, xmax, degree=5):
  s = s[:int(math.ceil(xmax/dx))]
  x = np.arange(0,xmax,dx)
  x += dx/2
  x = x[:len(s)]
  p = np.polyfit(x, s, degree)
  return np.polyval(p, 0.0)

full_volume_data = np.loadtxt(fr'{dir}/' + file + '.vol',skiprows=1)

sdf_data_init = np.loadtxt(fr'{dir}/' + file + '-init.sdf')
sdf_data_eq = np.loadtxt(fr'{dir}/' + file + '-eq.sdf')
sdf_data_npt = np.loadtxt(fr'{dir}/' + file + '-npt.sdf')

sdf_t_init = sdf_data_init[:,0]
sdf_data_init = sdf_data_init[:,1:]

sdf_t_eq = sdf_data_eq[:,0]
sdf_data_eq = sdf_data_eq[:,1:]

sdf_t_npt = sdf_data_npt[:,0]
sdf_data_npt = sdf_data_npt[:,1:]

volume_data_init = []
t_steps = list(map(int,sdf_t_init))
for vol in full_volume_data:
    if int(vol[0]) in t_steps:
        volume_data_init.append(float(vol[1]))
volume_data_init = np.array(volume_data_init)
N = full_volume_data[0,2]
n = len(sdf_data_init)
dim = 2
xmax = 0.001
dx = 1e-7

s0 = []
for i in range(n):
    s0.append(extrapolate(sdf_data_init[i], dx, xmax, degree = 5))
s0 = np.array(s0)
N = len(system.particles)
n_density = N/volume_data_init
betaP = n_density*(1.0 + s0/(2*dim))

with open(fr'{dir}/' + file + "-init.pressure",'w') as pv:
    pv.writelines('timestep:    volume:     betaP:\n')
    for i in range(len(t_steps)):
        pv.writelines(str(t_steps[i]) + ' ' + str(volume_data_init[i]) + '  ' + str(betaP[i]) + '\n')

plt.plot(t_steps,betaP,'--')
plt.title('Initial')
plt.savefig(fr"{dir}/Initial")
plt.show()

volume_data_eq = []
t_steps = list(map(int,sdf_t_eq))
for vol in full_volume_data:
    if int(vol[0]) in t_steps:
        volume_data_eq.append(float(vol[1]))
volume_data_eq = np.array(volume_data_eq)
N = full_volume_data[0,2]
n = len(sdf_data_eq)
dim = 2
xmax = 0.001
dx = 1e-7

s0 = []
for i in range(n):
    s0.append(extrapolate(sdf_data_eq[i], dx, xmax, degree = 5))
s0 = np.array(s0)
N = len(system.particles)
n_density = N/volume_data_eq
betaP = n_density*(1.0 + s0/(2*dim))

with open(fr'{dir}/' + file + "-eq.pressure",'w') as pv:
    pv.writelines('timestep:    volume:     betaP:\n')
    for i in range(len(t_steps)):
        pv.writelines(str(t_steps[i]) + ' ' + str(volume_data_eq[i]) + '  ' + str(betaP[i]) + '\n')

plt.plot(t_steps,betaP,'--')
plt.title('Equilibriate')
plt.savefig(fr"{dir}/Equilibriate")
plt.show()

volume_data_npt = []
t_steps = list(map(int,sdf_t_npt))
for vol in full_volume_data:
    if int(vol[0]) in t_steps:
        volume_data_npt.append(float(vol[1]))
volume_data_npt = np.array(volume_data_npt)
N = full_volume_data[0,2]
n = len(sdf_data_npt)
dim = 2
xmax = 0.005
dx = 1e-6

s0 = []
for i in range(n):
    s0.append(extrapolate(sdf_data_npt[i], dx, xmax, degree = 5))
s0 = np.array(s0)
N = len(system.particles)
n_density = N/volume_data_npt
betaP = n_density*(1.0 + s0/(2*dim))

with open(fr'{dir}/' + file + "-npt.pressure",'w') as pv:
    pv.writelines('timestep:    volume:     betaP:\n')
    for i in range(len(t_steps)):
        pv.writelines(str(t_steps[i]) + ' ' + str(volume_data_npt[i]) + '  ' + str(betaP[i]) + '\n')

plt.plot(N * tot_vol/volume_data_npt,betaP,'--')
plt.title('NPT - betaP')
plt.savefig(fr"{dir}/NPT")
plt.show()

os.system(fr'mv {file}.log {dir}')

