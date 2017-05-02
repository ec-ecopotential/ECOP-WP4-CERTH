# This python despeckling  module was implemented by CERTH within the framework of the ECOPOTENTIAL project.
# This program is free software. It comes without any warranty, to the extent permitted by applicable law. 
# You can redistribute it and/or modify it freely.


import gdal;
import cv2;
import numpy;
import math;
import sys;
import getopt;
import osr;
import gdalconst;
import os.path;
import array;

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
	
	
	
def main(argv):	
		input_guided_image = '';
		input_guidance_image ='';
		input_cloud_maskl = '';
		
		try:
			opts, args = getopt.getopt(argv,"h:i:g:m:r:e:o:",["ifile=","gfile=","radius=","smoothness_parameter=","mfile=","ofile="])
		except getopt.GetoptError:
			print 'Python1.py -i <input_file> -g <guidancefile> -r <radius> -e <smoothness_parameter> -m <cloud_mask> -o <outputfile>'
			try:
				opts, args = getopt.getopt(argv,"h i:g:r:e:o:",["ifile=","gfile=","radius=","smoothness_parameter=","ofile="])
			except getopt.GetoptError:
				print ' EITHER: Python1.py -i <input_file> -g <guidancefile> -o <outputfile>'
				sys.exit(2)					
		for opt, arg in opts:
			if opt == '-h':
				print 'Guided_filter_despeckling.py -i <inputfile> -g <guidancefile> -r <radius> -e <smoothness_parameter> -m <cloud_mask> -o <outputfile>'
				sys.exit()
			elif opt in ("-i", "--ifile"):
				inputfile=arg
			elif opt in ("-g", "--gfile"):
				guidancefile=arg
			elif opt in ("-r", "--radius"):
				r1=arg	
			elif opt in ("-e", "--smoothness_parameter"):
				eps1=arg		
			elif opt in ("-m", "--mfile"):
				maskfile=arg
			elif opt in ("-o","--ofile"):
				outputfile=arg
  
		if 'inputfile' not in locals(): #Check if the source file of the guided image was provided as input
			print "Source file of guided image was not provided."
			sys.exit(7)
		if 'guidancefile' not in locals(): #Check if the source file of the guidance image was provided as input
			print "Source file of guidance image was not provided."
			sys.exit(7)
		if 'outputfile' not in locals(): #Check if the filename of the output despeckled image was provided as input
			print "Output filename was not provided."
			sys.exit(7)
		if 'r1' not in locals(): #Set default radious value to 5
			r1=3;	
		if 'eps1' not in locals(): #Set default smoothness parameter value to 0.01
			eps1=0.001
				
		if os.path.isfile(guidancefile)==0: #Check the existence of the guidance source file and return a message when it does not exist
			print "Guidance image source file does not exist, please provide a new filename"
			sys.exit(7)
			
			
		print "Despeckling Started."
		
		Guidance_Image=gdal.Open(guidancefile) # Load guidance image
		
		if Guidance_Image.RasterCount !=3:
			print "Guidance image source file does not have three bands, please provide a new filename"
			sys.exit(7)
		
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
				
		if os.path.isfile(inputfile)==0: #Check the existence of the guidance source file and return a message when it does not exist
			print "Guided source file does not exist"
			sys.exit(7)
		
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
		   print "Projection systems between guided image and guidance image are not the same."
		   sys.exit(3)	
		   
		#print geoTrans_guided.IsSameGeogCS(geoTrans_guidance)			
		#kernel = numpy.ones((5,5),numpy.float32)/25
		#dst = cv2.filter2D(data,-1,kernel);
	
		r = int(r1);
		#print "Convert string to int %d" %r
		eps = float(eps1);
		#print "Convert string to float %f" %eps
		#eps = 0.01;
		
		#small = cv2.resize(Median_filter, (0,0), fx=0.10, fy=0.10) 
		#cv2.imshow('dst_rt', small)
		#cv2.waitKey(0)
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
			
		if 'maskfile' in locals(): #Check if the source file of the guidance image was provided as input
			if os.path.isfile(maskfile)==0: #Check the existence of the guidance source file and return a message when it does not exist
				print "Paparidis2"
				print "Cloud mask source file does not exist, cloud masking will not be considered in despeckling"							
				sys.exit(7)

			else:			    
				MASK_Cloud=gdal.Open(maskfile) # Load guidance image
				wkt_cloud_mask = MASK_Cloud.GetProjection() # Retrieve projection system of guided image into well known text (WKT) format and save it in wkt_guided
				srs_cloud_mask = osr.SpatialReference(wkt_guided) # Convert WTK format to OpenGIS Spatial Reference System object and save it in srs_guided
				if srs_cloud_mask.IsSame(srs_guided)==0:	#compare srs_cloud_mask with srs_guided in order to check if cloud mask and guided images have the same projection system
					print "Projection systems between cloud mask image and guided image are not the same."
					sys.exit(3)	
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
						print "Cloud mask file is not binary"
						sys.exit(3)
				
		#	where_are_NaNs = numpy.isnan(data_MASK)		
		#	data_MASK[where_are_NaNs]=1
		#	where_are_zeros = data_MASK==0	
		#	sum_of_zero_points=len(where_are_zeros)
		#	where_are_ones = data_MASK==1
		#	sum_of_one_points=len(where_are_ones)
		#	total_points=sum_of_one_points+sum_of_zero_points
		
		#	print "Gio is a legend %d" %total_points
		
		#small0 = cv2.resize(Ir, (0,0), fx=0.25, fy=0.25) 
		#cv2.imshow('Despeckle', small0)			
		#small1 = cv2.resize(q/255, (0,0), fx=0.25, fy=0.25) 
		#cv2.imshow('Despeckled Image', small1)
		#small2 = cv2.resize(Median_filter, (0,0), fx=0.25, fy=0.25) 
		#cv2.imshow('Median Image', small2)
		#cv2.waitKey(0)
		#print "Kia %f" % numpy.amin(q)
		#print "Kia1 %f" % (numpy.amax(q)*r)
		#print "Kia %f" % numpy.amin(p)
		#print "Kia1 %f" % numpy.amax(p)
		#print "eps1 %f" % eps
		q=(maximum_SAR-minimum_SAR)*q-minimum_SAR;
		
		format= 'GTiff'
		driver= gdal.GetDriverByName(format) #Generate an object of type Geotiff
		dst_ds = driver.Create(outputfile, Guided_Image.RasterXSize, Guided_Image.RasterYSize, 1, gdal.GDT_Float32) #Create a raster of type Geotiff with dimension Guided_Image.RasterXSize x Guided_Image.RasterYSize, with one band and datatype of GDT_Float32
		if dst_ds is None: #Check if outputfile can be saved
			print 'Could not save output file %s, path does not exist.' % outputfile
			sys.exit(4)
			
		dst_ds.SetGeoTransform(geoTrans_guided) # Set the Geoinformation of the output file the same as the one of the guidance image
		dst_ds.SetProjection (wkt_guidance)
		dst_ds.GetRasterBand(1).WriteArray(q) # Save the raster into the output file	
		dst_ds.FlushCache()  # Write to disk.			
		#wkt_output= dst_ds.GetProjection();
		
		print "Despeckling Finished."
			
		#print wkt_guided;
		#print "Keno"
		#print wkt_output;
		#cv2.imwrite(outputfile,255*q);
		#cv2.imwrite(outputfile,255*p);
		#Saved_Image=gdal.Open("D:/deleteme.tif")
		#wkt_saved_Image = Saved_Image.GetProjection();
		#print wkt_saved_Image=    
	
if __name__ == "__main__":
	main(sys.argv[1:])	
	
else:
    print("Input arguments were not provided")


sys.exit();


	
