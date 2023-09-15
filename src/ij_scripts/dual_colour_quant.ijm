macro "extract_counting_quantification"{
ori=getTitle();
getDimensions(dummy,dummy,ch, slices, frames);
file8 = getDirectory("quantifs");
ROIsize=getNumber("With of ROI around the spot [px]",15);
//IJ.renameResults("All Spots statistics","Results");
r=ROIsize/2;

nR=nResults(); //number of rows in the Results table

f1=File.open(file8 +'f1.txt');
File.close(f1);
f2=File.open(file8 +'f2.txt');
File.close(f2);

//retreive track IDs from results table
print("\\Clear");
for(i=0;i<nR;i++){
			Z=getResult("Slice", i);
			X=getResult("X", i);
			Y=getResult("Y", i);
			selectWindow(ori);
			Stack.setSlice(Z);
			makeRectangle(X-r, Y-r, ROIsize, ROIsize);
			y=i+1;
			ori=getTitle();
			zb=Z-7;
			ze=Z+8;
			run("Duplicate...", "duplicate slices=zb-ze");
			saveAs("Tiff",file8 + ori +'crop'+ -y);
			run("Z Project...", " projection=[Max Intensity]");
			saveAs("Tiff",file8 + ori +'Max'+ -y); 	
			run("Select All");
			profile= getProfile();
			Xvalues=Array.getSequence(lengthOf(profile));
			//print(lengthOf(profile));
			Array.getStatistics(profile,min,max,mean,std);
			
			InitGuesses=newArray(min,(profile[profile.length-1]-profile [0])/(lengthOf(profile)), max-min, round(profile.length)/2, 1);
			Fit.doFit( "y= a+b*x+c*exp(-(x-d)*(x-d)/e) ", Xvalues, profile, InitGuesses);
			Fit.plot;
			saveAs("Tiff", file8+ i+1 +"-PLot1"); 
			Sum1=0;
			Corrprofile=Array.copy(profile);
			for(j=0;j<lengthOf(profile);j++){
			Corrprofile[j]= profile[j]- Fit.p(0)-j*Fit.p(1);
			Sum1+=Corrprofile[j];
			}
			File.append(d2s(Sum1,2), file8 +'f1.txt');
			File.close(f1);
			close();
			
			run("Next Slice [>]");
			run("Select All");
			profile= getProfile();
			Xvalues=Array.getSequence(lengthOf(profile));
			//print(lengthOf(profile));
			Array.getStatistics(profile,min,max,mean,std);
			
			InitGuesses=newArray(min,(profile[profile.length-1]-profile [0])/(lengthOf(profile)), max-min, round(profile.length)/2, 1);
			Fit.doFit( "y= a+b*x+c*exp(-(x-d)*(x-d)/e) ", Xvalues, profile, InitGuesses);
			Fit.plot;
			saveAs("Tiff", file8+ i+1 +"-PLot2"); 
			Sum2=0;
			Corrprofile=Array.copy(profile);
			for(j=0;j<lengthOf(profile);j++){
			Corrprofile[j]= profile[j]- Fit.p(0)-j*Fit.p(1);
			Sum2+=Corrprofile[j];
			}
			
			File.append(d2s(Sum2,2), file8 +'f2.txt');
			File.close(f2);
			close();
			saveAs("Jpeg",file8 + ori +'Max'+ -y);
			close();
			run("Split Channels");
			run("Z Project...", " projection=[Max Intensity]");
			setMinAndMax(100, 300);
			saveAs("Tiff",file8 + ori + -y + "-C1");
			saveAs("Jpeg",file8 + ori + -y + "-C1");
			close();
			close();
			run("Z Project...", " projection=[Max Intensity]");
			setMinAndMax(100, 300);
			saveAs("Tiff",file8 + ori + -y + "-C2");
			saveAs("Jpeg",file8 + ori + -y + "-C2");
			close();
			close();
		}	
run("From ROI Manager");
saveAs("TIFF", file8+ori+"-spots");

}
