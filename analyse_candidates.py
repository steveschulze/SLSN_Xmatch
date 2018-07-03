import	argparse
from 	astropy import table
from	astropy import coordinates as coord
from	astropy import units as u
from 	astroquery import vizier
import	json
from 	misc import bcolors
import 	multiprocessing as mp
import	numpy as np 
import 	warnings

warnings.filterwarnings("ignore")

### Routines

def wrapper(x):
	return check_object(x, catalog_prop)

def vizier_query(OBJECT, PHOTCAT, CATPROP, ROW_LIMIT=-1, RADIUS=10. * u.arcmin):

	v 		= vizier.Vizier(columns = CATPROP[PHOTCAT]['INPUT'], catalog = CATPROP[PHOTCAT]['CATID'], row_limit = ROW_LIMIT)
	result	= v.query_region(coord.SkyCoord(OBJECT['ra'], OBJECT['dec'], unit=u.deg), radius = RADIUS)
	if len(result) > 0:
		try:
			result	= result[CATPROP[PHOTCAT]['CATID_OUT']][CATPROP[PHOTCAT]['OUTPUT']]
			return result
		except:
			return table.Table(names=CATPROP[PHOTCAT]['OUTPUT'])
	else:
		return table.Table(names=CATPROP[PHOTCAT]['OUTPUT'])

def check_object(OBJECT, CATALOGS):

	if args.verbose:
		print(bcolors.HEADER  + bcolors.BOLD + OBJECT['name'] + bcolors.ENDC)
		print(bcolors.WARNING + 'RA, DEC = {ra:.6f}, {dec:.6f}'.format(ra=OBJECT['ra'], dec=OBJECT['dec']) + bcolors.ENDC)
		print('\n')

	### Query SDSS

	results_sdss		= vizier_query(OBJECT, 'SDSS_DR12', CATALOGS, RADIUS=2. * u.arcsec)

	if args.verbose:
		print(bcolors.OKBLUE + 'Unfiltered SDSS output' + bcolors.ENDC)
		print (results_sdss)
		print('\n')

	flag_qso= False
	spec_z	= np.array([np.nan, np.nan])
	photo_z	= np.array([np.nan, np.nan])

	if len(results_sdss) > 0:
		results_sdss.sort('_r')
		results_sdss				= results_sdss[(results_sdss['mode'] == 1)]# & (results_sdss['q_mode'] == '+')]

		if args.verbose:
			print(bcolors.OKGREEN + 'Rejecting bad data' + bcolors.ENDC)
			print(results_sdss)
			print('\n')

		if len(results_sdss)		> 0:
	
			results_sdss['FLAG_QSO']	= [True if 'AGN' in x or 'BROADLINE' in x else False for x in results_sdss['subCl']]
			
			if args.verbose:
				print(bcolors.OKGREEN + 'Final SDSS output' + bcolors.ENDC)
				print(results_sdss)
				print('\n')

			if not results_sdss['FLAG_QSO'][0]:
				spec_z 	= results_sdss[['zsp', 'e_zsp']][0] if results_sdss['zsp'] > 0 else np.array([np.nan, np.nan])
				photo_z	= results_sdss[['zph', 'e_zph']][0] if results_sdss['zph'] > 0 else np.array([np.nan, np.nan])
				
			else:
				flag_qso= True

	### Query GAIA DR2

	results_gaia		= vizier_query(OBJECT, 'GAIA_DR2', CATALOGS, RADIUS=40. * u.arcsec)

	if args.verbose:
		print(bcolors.OKBLUE + 'Unfiltered Gaia DR2 output' + bcolors.ENDC)
		print (results_gaia)
		print('\n')

	flag_star			= False

	if len(results_gaia) > 0:

		results_gaia	= results_gaia[(results_gaia['Plx'] > 0) | (results_gaia['pmRA'] > 0) | (results_gaia['pmDE'] > 0)]

		# print (results_gaia)

		if len(results_gaia) > 0:

			results_gaia.sort('_r')

			# Identify stars
			# Thresholds for parallax and proper-motion measurements

			sig_plx			= 3
			sig_pmra		= 3
			sig_pmdec		= 3

			# If the significance of either the parallax or one of the proper-motion measurements exceeds N sigma, 

			mask_star		= [True if x['Plx'] / x['e_Plx'] >= sig_plx or x['pmRA'] / x['e_pmRA'] >= sig_pmra or x['pmRA'] / x['e_pmRA'] >= sig_pmra else False for x in results_gaia]
			results_gaia	= results_gaia[mask_star]

			if args.verbose:
				print(bcolors.OKGREEN + 'After removing all objects with either insignificant parallax or proper-motion measurements' + bcolors.ENDC)
				print (results_gaia)
				print('\n')

			# Remove objects that are too close to stars
			# Ensuring good observability and avoid artefacts from saturated stars

			if len(results_gaia) > 0:

				results_gaia['DISTANCE_NORM']	= 1.8 + 0.6 * np.exp( (20 - results_gaia['Gmag']) / 2.05)
				results_gaia['FLAG_PROX']		= [True if x['DISTANCE_NORM'] > x['_r'] and 11 <= x['Gmag'] <= 19 else False for x in results_gaia]

				if args.verbose:
					print(bcolors.OKGREEN + 'Flag stars that are too close to the transient position for their given brightness.' + bcolors.ENDC)
					print (results_gaia)
					print('\n')

				if any(results_gaia['FLAG_PROX'] == True):
					flag_star	= True

	### Query Milliquas

	results_milliquas		= vizier_query(OBJECT, 'Milliquas', CATALOGS, RADIUS=2. * u.arcsec)

	if args.verbose:
		print(bcolors.OKBLUE + 'Unfiltered Milliquas output' + bcolors.ENDC)
		print (results_milliquas)
		print('\n')

	if flag_star or flag_qso:
		return None
	else:
		return {'NAME': OBJECT['name'], 'SPEC_Z': spec_z, 'PHOTO_Z': photo_z, 'MILLIQUAS': results_milliquas}

### Input parameters

parser				 = argparse.ArgumentParser(description='Analyse output from the Marshal Scanning page')

# Object properties

parser.add_argument('--json',			type	= str,
										help	= 'JSON input file',
										required= True)

parser.add_argument('--ncpu',			type	= int,
										help	= 'Number of CPUs (Default: Core x Threading - 1)',
										default	= mp.cpu_count()-1)

parser.add_argument('--verbose',		action	= 'store_true',
										default	= False,
										help	= 'Verbose')

args 				= parser.parse_args()

# Read JSON file
with open(args.json) as data_file:
	candidates = json.load(data_file)

# Catalogue properties

catalog_prop							= {}
catalog_prop['GAIA_DR2']				= {}
catalog_prop['GAIA_DR2']['INPUT'] 		= ['_r', '*']
catalog_prop['GAIA_DR2']['OUTPUT'] 		= ['_r', 'RA_ICRS', 'DE_ICRS', 'Gmag', 'Plx', 'e_Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE']
catalog_prop['GAIA_DR2']['CATID']		= "I/345"
catalog_prop['GAIA_DR2']['CATID_OUT']	= catalog_prop['GAIA_DR2']['CATID'] + '/gaia2'

catalog_prop['SDSS_DR12']				= {}
catalog_prop['SDSS_DR12']['INPUT']		= ['*', '_r', 'subCl', 'e_zsp', 'e_zph']
catalog_prop['SDSS_DR12']['OUTPUT'] 	= ['_r', 'RA_ICRS', 'DE_ICRS', 'mode', 'q_mode', 'class', 'subCl', 'zsp', 'zph', 'e_zsp', 'e_zph']
catalog_prop['SDSS_DR12']['CATID']		= "V/147"
catalog_prop['SDSS_DR12']['CATID_OUT']	= catalog_prop['SDSS_DR12']['CATID'] + '/sdss12'

catalog_prop['Milliquas']				= {}
catalog_prop['Milliquas']['INPUT']		= ['*', '_r']
catalog_prop['Milliquas']['OUTPUT'] 	= ['_r', 'RAJ2000', 'DEJ2000', 'Name', 'Cl', 'Qpct', 'z']
catalog_prop['Milliquas']['CATID']		= "VII/280"
catalog_prop['Milliquas']['CATID_OUT']	= catalog_prop['Milliquas']['CATID'] + '/catalog'

# Catalog queries

if args.ncpu		> 1:
	pool 		= mp.Pool(processes=7)
	results 	= pool.map(wrapper, candidates)
else:
	results	= [check_object(x, catalog_prop) for x in candidates]

# Process output

results	= [x for x in results if x is not None]

# Text output

print('\n')
print('Number of candidates: {number:g}'.format(number=len(candidates)))
print('Number of filtered candidates: {number:g}'.format(number=len(results)))
print('\n')

for result in results:
	print(bcolors.HEADER + 'Object:  {name}'.format(name=result['NAME']) + bcolors.ENDC)
	print('spec-z:  {value} +/- {error}'.format(value=result['SPEC_Z'][0], error=result['SPEC_Z'][1]))
	print('photo-z: {value} +/- {error}'.format(value=result['PHOTO_Z'][0], error=result['PHOTO_Z'][1]))
	print('Milliquas output')
	print(result['MILLIQUAS'])
	print('')