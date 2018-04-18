# This python despeckling  module was implemented by CERTH within the framework of the ECOPOTENTIAL project.
# This program is free software. It comes without any warranty, to the extent permitted by applicable law. 
# You can redistribute it and/or modify it freely.


import gdal
import cv2
import numpy
import math
import osr
import os.path
import zipfile
import time

supported_img_formats = [".tif", ".tiff"]
main_directory = "data/"
zip_file = main_directory + 'despeckling_input.zip'
settings_filename = "despeckling_settings.txt"
#
identifiers = ["-i", "-g", "-m", "-r", "-e" , "-o"]
idsDict = dict.fromkeys(identifiers)

def range_0_1(data): #Values are clamped to the 1st and 99th percentile and scaled to 0-1
	where_are_NaNs = numpy.isnan(data)
	data[where_are_NaNs]=0;
	minimum_percentile=numpy.percentile(-data,99)
	#print "minimum_percentile %f " % minimum_percentile;
	maximum_percentile=numpy.percentile(data,99)
	#print "maximum_percentile %f " % maximum_percentile
	data[data < minimum_percentile] = minimum_percentile;
	data[data > maximum_percentile] = maximum_percentile;
	data1=data-minimum_percentile;
	data1=numpy.divide(data1,numpy.amax(data1));
	return data1
	
def range_1(data): #Values are clamped to the 1st and 99th percentile and scaled to 0-1
	
	minimum_percentile=numpy.amin(data)
	#print "minimum_percentile %f " % minimum_percentile;
	maximum_percentile=numpy.amax(data)
	#print "maximum_percentile %f " % maximum_percentile	
	data1=data-minimum_percentile;
	data1[data1 < 0]=0;
	data1=numpy.divide(data1,numpy.amax(data1));
	return data1, minimum_percentile, maximum_percentile
	

#compress the contents of a given directory to a zip file
def zip(src, dst):
    zf = zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED)
    abs_src = os.path.abspath(src)
    for dirname, subdirs, files in os.walk(src):
        for filename in files:
            absname = os.path.abspath(os.path.join(dirname, filename))
            arcname = absname[len(abs_src) + 1:]
            zf.write(absname, arcname)
    zf.close()


#read text file and validates that referenced files from fragmentation_settings.txt exist in the archive
def readZipAttributes(zip_file_path):
	with zipfile.ZipFile(zip_file_path) as archive_input:
		#check if zip file is corrupt
		try:
			ret = archive_input.testzip()
			if ret is not None:
				print ("First bad file in zip: %s" % ret)
				writeReport("First bad file in zip: %s" % ret)
				return -1
			else:
				print "Zip file is good."
		except:
			print ("Zip is corrupted1!")
			writeReport("Zip is corrupted!")
			return -1
		if settings_filename in archive_input.namelist():
			print("Configuration file '" + settings_filename + "' detected.")
			txtfile = archive_input.read(settings_filename)
			parameters_list = txtfile.split()

			#parameters = " ".join(txtfile.split())
			for identifier in identifiers:
				indices_found = [i for i, x in enumerate(parameters_list) if x == identifier]
				if len(indices_found) > 1:
					print("Identifier '" + identifier +"' appears multiple times")
					return -1
				elif len(indices_found) == 1:
					try:
					    idsDict[identifier] = parameters_list[indices_found[0] + 1]
					except IndexError:
						print("Given Identifier '" + identifier +"' is followed by no value")
						return -1				

			if idsDict["-i"] is not None and idsDict["-i"] not in archive_input.namelist():
				print ("SAR raster file %s does not exist within the archive, although declared in %s" % (idsDict["-i"], settings_filename))
				return -1
			if idsDict["-g"] is not None and idsDict["-g"] not in archive_input.namelist():
				print ("RGB raster file %s does not exist within the archive, although declared in %s" % (idsDict["-g"], settings_filename))
				return -1
			if idsDict["-m"] is not None and idsDict["-m"] not in archive_input.namelist():
				print ("Cloud Mask raster file %s does not exist within the archive, although declared in %s" % (idsDict["-m"], settings_filename))
				return -1	
		return idsDict
	
#writes a report with a given string (message) to a fixed output file
def writeReport(single_msg):
	with open(main_directory + 'outputs/io_report'+'.txt','w') as f:
		f.write(single_msg)
	
def main(dict_args):
	
	if not os.path.exists(main_directory + 'outputs'):
		os.mkdir(main_directory + 'outputs')

	start_time = time.time()
					
	#assign variables
	inputfile = dict_args["-i"]
	guidancefile = dict_args["-g"]
	r1 = dict_args["-r"]
	eps1 = dict_args["-e"]
	maskfile = dict_args["-m"]
	outputfile = dict_args["-o"]
	

	#validate fragmentation_settings.txt
	if inputfile is None: #Check if the source file of the guided image was provided as input
		print("Source file of guided image was not provided.")
		return
	if guidancefile is None: #Check if the source file of the guidance image was provided as input
		print("Source file of guidance image was not provided.")
	if outputfile is None: #Check if the filename of the output despeckled image was provided as input
		print("Output filename was not provided.")
	if r1 is None: #Set default radious value to 3
		r1 = 3
		print("Radius was not provided.")		
	if eps1 is None: #Set default smoothness parameter value to 0.001
		eps1 = 0.001
		print("Smoothness parameter was not provided.")
	
	print("Radius = %s and Smoothness = %s" % (r1, eps1))


	#check if given rasters are in supported format		
	input_file_filename_only, input_file_file_extension = os.path.splitext(inputfile)
	guidancefile_filename_only, guidancefile_file_extension = os.path.splitext(guidancefile)
	#maskfile could be None (it's optional)
	maskfile_filename_only, maskfile_file_extension = os.path.splitext(maskfile)


	if input_file_file_extension not in supported_img_formats or guidancefile_file_extension not in supported_img_formats or (maskfile is not None and maskfile_file_extension not in supported_img_formats):
		print ("Critical error! Provided rasters's extenion is not supported. Valid extensions are: %s" % supported_img_formats)
		writeReport("Critical error! Provided rasters's extenion is not supported. Valid extensions are: %s" % supported_img_formats)		
		return -1
	else:
		print ("Formats are supported.")		
		
	print "Despeckling Started."
	
	guidancefile = '/vsizip/' + zip_file + '/%s' % guidancefile
	Guidance_Image=gdal.Open(guidancefile) # Load guidance image
	
	if Guidance_Image.RasterCount !=3:
		print ("Guidance image source file does not have three bands, please provide a new filename")
		writeReport("Guidance image source file does not have three bands, please provide a new filename")
		return -1
	
	srcband = Guidance_Image.GetRasterBand(1)
	data_R = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.float) # Save Red channel to data_R
	srcband = Guidance_Image.GetRasterBand(2)
	data_G = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.float) # Save Green channel to data_G
	srcband = Guidance_Image.GetRasterBand(3)
	data_B = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.float) # Save Blue channel to data_B	
	minimum_percentile=numpy.percentile(-data_R,99)
	median_value=numpy.median(data_R);
	difference = numpy.abs(minimum_percentile-median_value)		
			
	if difference>median_value*100:
		data_R = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.integer) # Save Red channel to data_R
		srcband = Guidance_Image.GetRasterBand(2)
		data_G = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.integer) # Save Green channel to data_G
		srcband = Guidance_Image.GetRasterBand(3)
		data_B = srcband.ReadAsArray(0, 0, Guidance_Image.RasterXSize, Guidance_Image.RasterYSize).astype(numpy.integer) # Save Blue channel to data_B				
	
	geoTrans_guidance = Guidance_Image.GetGeoTransform() # Retrieve Geoinformation of guidance image and save it in geoTrans_guidance
	wkt_guidance = Guidance_Image.GetProjection() # Retrieve projection system of guidance image into well known text (WKT) format and save it in wkt_guidance
	srs_guidance=osr.SpatialReference(wkt_guidance) # Convert WTK format to OpenGIS Spatial Reference System object and save it in srs_guidance

	data_R1=range_0_1(data_R) #Scaled Red channel is saved in data_R1
	data_G1=range_0_1(data_G) #Scaled Green channel is saved in data_G1
	data_B1=range_0_1(data_B) #Scaled Blue channel is saved in data_B1
	
	print "RGB Guidance image has been loaded."
			
	inputfile = '/vsizip/' + zip_file + '/%s' % inputfile
	
	Guided_Image=gdal.Open(inputfile) #Load guided image
	srcband = Guided_Image.GetRasterBand(1);
	data_S = srcband.ReadAsArray(0, 0, Guided_Image.RasterXSize, Guided_Image.RasterYSize).astype(numpy.float) #Save guided image to data_S
	minimum_percentile=numpy.percentile(-data_S,99)
	median_value=numpy.median(data_S);
	difference = numpy.abs(minimum_percentile-median_value)
	
	#print "Median_value %f and difference %f" % (median_value, difference)
	if difference>median_value*100:
		data_S = srcband.ReadAsArray(0, 0, Guided_Image.RasterXSize, Guided_Image.RasterYSize).astype(numpy.integer) 
	
	data_SAR,minimum_SAR,maximum_SAR = range_1(data_S) #Scaled guided image is stored in data_SAR
	
	print "SAR image to be despeckled has been loaded."
			
		
	geoTrans_guided = Guided_Image.GetGeoTransform() # Retrieve Geoinformation of guided image and save it in geoTrans_guidance
	wkt_guided = Guided_Image.GetProjection() # Retrieve projection system of guided image into well known text (WKT) format and save it in wkt_guided	
	srs_guided=osr.SpatialReference(wkt_guided) # Convert WTK format to OpenGIS Spatial Reference System object	and save it in srs_guided


	#print srs_guided
	#print srs_guidance
	
	if srs_guided.IsSame(srs_guidance)==0:	#compare srs_guidance with srs_guided in order to check if guided and guidance images have the same projection system
		print("Projection systems between guided image and guidance image are not the same.")
		writeReport("Projection systems between guided image and guidance image are not the same.")
		return -1	

	r = int(r1);
	#print "Convert string to int %d" %r
	eps = float(eps1);
	#print "Convert string to float %f" %eps

	
	print  "Guided filtering despeckling process has initiated."
    
	#****************** Start of applying the guided filter approach ******************************
	r2=2*r+1;
	
	Ir=data_R1
	del data_R1
	Ig=data_G1;
	del data_G1;
	Ib=data_B1;
	del data_B1;
	p=data_SAR;
	del data_SAR;

	Ir_mean = cv2.blur(Ir,(r2, r2));
		
	
	Ig_mean = cv2.blur(Ig,(r2, r2));
	Ib_mean = cv2.blur(Ib,(r2, r2));

	p_mean = cv2.blur(p,(r2, r2));

	Ipr_mean = cv2.blur(Ir * p, (r2, r2));
	Ipg_mean = cv2.blur(Ig * p, (r2, r2));
	Ipb_mean = cv2.blur(Ib * p, (r2, r2));

	Ipr_cov = Ipr_mean - Ir_mean * p_mean;
	Ipg_cov = Ipg_mean - Ig_mean * p_mean;
	Ipb_cov = Ipb_mean - Ib_mean * p_mean;
	
	del Ipr_mean
	del Ipg_mean
	del Ipb_mean

	#Sigma + eps*eye(3)
	Irr_var = cv2.blur(Ir * Ir, (r2, r2)) - Ir_mean * Ir_mean + eps; 
	Irg_var = cv2.blur(Ir * Ig, (r2, r2)) - Ir_mean * Ig_mean;
	Irb_var = cv2.blur(Ir * Ib, (r2, r2)) - Ir_mean * Ib_mean;
	Igg_var = cv2.blur(Ig * Ig, (r2, r2)) - Ig_mean * Ig_mean + eps;
	Igb_var = cv2.blur(Ig * Ib, (r2, r2)) - Ig_mean * Ib_mean;
	Ibb_var = cv2.blur(Ib * Ib, (r2, r2)) - Ib_mean * Ib_mean + eps;
			

	Irr_inv = Igg_var * Ibb_var - Igb_var * Igb_var;
	Irg_inv = Igb_var * Irb_var - Irg_var * Ibb_var;
	Irb_inv = Irg_var * Igb_var - Igg_var * Irb_var;
	Igg_inv = Irr_var * Ibb_var - Irb_var * Irb_var;
	Igb_inv = Irb_var * Irg_var - Irr_var * Igb_var;
	Ibb_inv = Irr_var * Igg_var - Irg_var * Irg_var;
   

	I_cov = Irr_inv * Irr_var + Irg_inv * Irg_var + Irb_inv * Irb_var;

	Irr_inv /= I_cov;
	Irg_inv /= I_cov;
	Irb_inv /= I_cov;
	Igg_inv /= I_cov;
	Igb_inv /= I_cov;
	Ibb_inv /= I_cov;

	del Irr_var
	del Irg_var
	del Irb_var
	del Igg_var
	del Igb_var
	del Ibb_var

	ar = Irr_inv * Ipr_cov + Irg_inv * Ipg_cov + Irb_inv * Ipb_cov;
	ag = Irg_inv * Ipr_cov + Igg_inv * Ipg_cov + Igb_inv * Ipb_cov;
	ab = Irb_inv * Ipr_cov + Igb_inv * Ipg_cov + Ibb_inv * Ipb_cov;
	b = p_mean - ar * Ir_mean - ag * Ig_mean - ab * Ib_mean;
	
	del Irr_inv
	del Irg_inv
	del Irb_inv
	del Igg_inv
	del Igb_inv
	del Ibb_inv

	ar_mean = cv2.blur(ar, (r2, r2));
	ag_mean = cv2.blur(ag, (r2, r2));
	ab_mean = cv2.blur(ab, (r2, r2));
	b_mean = cv2.blur(b, (r2, r2));

	q = (ar_mean * Ir + ag_mean * Ig + ab_mean * Ib + b_mean); # q is the despeckled guided image
	
	#******************* End of applying the guided filter approach ******************************	
	print "Guided filtering despeckling process has finished."
		
	if maskfile is not None: #Check if the source file of the guidance image was provided as input	
		maskfile = '/vsizip/' + zip_file + '/%s' % maskfile	    
		MASK_Cloud=gdal.Open(maskfile) # Load guidance image
		wkt_cloud_mask = MASK_Cloud.GetProjection() # Retrieve projection system of guided image into well known text (WKT) format and save it in wkt_guided
		srs_cloud_mask = osr.SpatialReference(wkt_guided) # Convert WTK format to OpenGIS Spatial Reference System object and save it in srs_guided
		if srs_cloud_mask.IsSame(srs_guided)==0:	#compare srs_cloud_mask with srs_guided in order to check if cloud mask and guided images have the same projection system
			print ("Projection systems between cloud mask image and guided image are not the same.")
			writeReport("Projection systems between cloud mask image and guided image are not the same.")
			return -1	
		else:
			srcband = MASK_Cloud.GetRasterBand(1)
			data_MASK = srcband.ReadAsArray(0, 0, MASK_Cloud.RasterXSize, MASK_Cloud.RasterYSize).astype(numpy.integer)
			counter=0;
			counter1=0;
			for x in range(0, MASK_Cloud.RasterXSize):
				for y in range(0, MASK_Cloud.RasterYSize):
					if data_MASK[y,x]==0:
						counter=counter+1
					elif data_MASK[y,x]==1:
						counter1=counter1+1
		
			if (counter+counter1) == (MASK_Cloud.RasterXSize*MASK_Cloud.RasterYSize):	
				print "Median filtering is performed for areas covered by clouds."
				kernel = numpy.ones((2*r+1,2*r+1),numpy.float32)/((2*r+1)*(2*r+1))
				Median_filter = cv2.filter2D(p,-1,kernel)										
				q[data_MASK==1]=Median_filter[data_MASK==1]
			else:
				print ("Cloud mask file is not binary")
				writeReport("Cloud mask file is not binary")
				return -1
			
	q=(maximum_SAR-minimum_SAR)*q-minimum_SAR;
	
	format= 'GTiff'
	driver= gdal.GetDriverByName(format) #Generate an object of type Geotiff
	dst_ds = driver.Create(main_directory + 'outputs/' + outputfile, Guided_Image.RasterXSize, Guided_Image.RasterYSize, 1, gdal.GDT_Float32) #Create a raster of type Geotiff with dimension Guided_Image.RasterXSize x Guided_Image.RasterYSize, with one band and datatype of GDT_Float32
	if dst_ds is None: #Check if outputfile can be saved
		print ('Could not save output file %s, path does not exist.' % outputfile)
		writeReport('Could not save output file %s, path does not exist.' % outputfile)
		return -1
		
	dst_ds.SetGeoTransform(geoTrans_guided) # Set the Geoinformation of the output file the same as the one of the guidance image
	dst_ds.SetProjection (wkt_guidance)
	dst_ds.GetRasterBand(1).WriteArray(q) # Save the raster into the output file	
	dst_ds.FlushCache()  # Write to disk.

	writeReport('Succesfully created %s' % outputfile)		

	elapsed_time = time.time() - start_time
	print elapsed_time
	
	print "Despeckling Finished."

if __name__ == "__main__":
	dict_params = readZipAttributes(zip_file)
	print(dict_params)
	if dict_params != -1:
		main(dict_params)

	# save result as zip
	print("Compressing results in a 'zip' file...")
	zip(main_directory + 'outputs/', main_directory + 'despeckling_output.zip')
	
	
else:
    print("Input arguments were not provided")



