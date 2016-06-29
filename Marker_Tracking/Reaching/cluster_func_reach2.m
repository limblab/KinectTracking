function [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func_reach2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 )
%This function clusters the points, and assigns them to markers.

%Inputs:
%t: the time point being looked at
%loc: The x,y,z coordinates of all the points (pixels of the color being looked at). 
%num_clust: The number of clusters to try to make (the number of total
%markers in the given section)
%prev_num_clust: The  number of markers detected in the previous frame
%dist_min: The minimum distance allowed between markers
%prev_meds: The locations of the markers in the previous frame
%medians: The locations of the markers at all previous times. Missing
%markers have NaN values.
%medians2: The locations of the markers at all previous times. Missing
%markers have the location at the last known time.

%Outputs:
%prev_num_clust: The number of markers detected in the current frame (which
%will be the number of markers detected in the previous frame, at the next
%time point).
%prev_meds: The locations of the markers detected in the current frame
%(which will be the previous locations of markers, at the next time point)
%medians/medians2: The locations of the markers at all times up through the current
%frame (the inputs are updated to include the current frame). The
%difference between medians and medians2 is the same as for the inputs.


%If there are no points (pixels of the color being looked at), set medians
%to have all NaNs. And set medians2 to be medians2 of the previous frame.
if isempty(loc)
    medians(:,:,t)=NaN;
    if t>1
        medians2(:,:,t)=medians2(:,:,t-1);
    else
        medians2(:,:,t)=0;%[0 0 0];
    end
else
    %We initiallly set the number of clusters to look for (curr_num_clust)
    %to be one greater than the true number of markers (num_clust). This is
    %because there are instances when some noise will get detected by the
    %extra cluster, and can then be removed. Without the extra cluster,
    %there are times when the noise is included as a cluster, and the true
    %marker is not included.
    curr_num_clust=num_clust+1;
    med_dist=zeros(1,num_clust); %Initialize med_dist so that  we will enter the while loop below
    
    if size(loc,1)<curr_num_clust %In case there are fewer points than clusters I'm fitting
            curr_num_clust=size(loc,1); %Set the number of clusters to be the number of points
    end
        
    %Each run within the while loop (below) clusters with "curr_num_clust" clusters.
    %At the end of the while loop, if the distance between any clusters is
    %too close (so that those clusters really should be combined), then we
    %reduce the number of clusters by 1, and rerun everything within the
    %while loop.
    %The conditions for doing another run through the while loop are:
    %1. If any of the distances between the markers (med_dist) are less
    %than dist_min OR
    %2. If any of the distances between the markers are NaNs (which happens if no points are assigned to a cluster).
    % Additionally, we don't run through the while loop if there if there
    % curr_num_clust=0 (if the previous run had curr_num_clust=1, meaning it used 1 cluster).
    while ((any(med_dist<dist_min) || any(isnan(med_dist))) && curr_num_clust>=1)      
        
        %1. Below, we will cluster the points. Depending on the situation, we
        %will sometimes use known starting points. Note that when we don't
        %have starting points, we use kmedoids rather than kmeans because
        %it is more robust.
        
        if curr_num_clust>num_clust %This is the first time through the loop
            %When curr_num_clust is 1 greater than prev_num_clust, we
            %cluster using the previous markers as starting points, with the
            %addition of an extra starting point at the first marker's location.
            if curr_num_clust==prev_num_clust+1 
                prev_meds2=[prev_meds(1,:); prev_meds];
                clust=kmeans(loc,curr_num_clust,'Start',prev_meds2);
            else
            %If curr_num_clust is >1 greater than prev_num_clust, we
            %cluster w/o starting points
                clust=kmedoids(loc,curr_num_clust);
            end
        else %When curr_num_clust <=num_clust 
        %If the curr_num_clust doesn't equal prev_num_clust cluster w/o
        %starting points
        if prev_num_clust~=curr_num_clust 
            clust=kmedoids(loc,curr_num_clust);
        else
        %If curr_num_clust equals prev_num_clust, we cluster w/ the
        %previous marker locations as starting points
            clust=kmeans(loc,curr_num_clust,'Start',prev_meds);
        end
        end
        
        %2. Record the median of each cluster (which is the marker
        %location)
    
        meds=zeros(curr_num_clust,3); %Initialize the medians (meds)
        for j=1:curr_num_clust %Loop through the clusters
            if nnz(clust==j)>1 %If there are more than 1 point in a cluster, find the median
                meds(j,:)=median(loc(clust==j,:));
            else if nnz(clust==j)==1 %If there is 1 point in a cluster, the median is that point's location
                    meds(j,:)=loc(clust==j,:);
                else %If there are no points in a cluster, set the median to NaNs.
                    meds(j,:)=NaN(1,3);
                end
            end
        end
        
        %3. Calculate the distance between markers
        med_dist=pdist(meds);
        
        %4. Lower the number of clusters by 1 for the next run in the while
        %loop if any of the distances between markers is too small or NaNs
        %(when there were no points in a cluster).
        if any(med_dist<dist_min) || any(isnan(med_dist))
            curr_num_clust=curr_num_clust-1;
        end
        
    end
    
    %5. Assign to markers (based on how close they are to previous marker locations - or could be some estimate of current locations)
    
    %If it's the first frame, just assign. Otherwise, match up with assignment
    %from previous frames. We match up (assign) by minimizing the sum squared
    %distance between the current markers and the markers from the previous
    %frame.
       
    if t==1 %First frame
        %For now, all of the clusters have to be present in the first frame analyzed
        medians(:,:,t)=0;%meds;
        medians2(:,:,t)=0;%meds;
        meds2=zeros(curr_num_clust,3);%meds;
        
    else %Not first frame
        temp=pdist2(medians2(:,:,t-1),meds); %Distance between all current markers and markers in previous frame
        temp2=temp.^2; %Distances squared
        
        %The above matrix of distances is a matrix of size num_clust x curr_num_clust
        % When num_clust=curr_num_clust, then we want to switch the columns
        % so that the values along the diagonals are minimal (and we just aim to minimize the trace of the "switched" matrix).
        %When num_clust~=curr_num_clust, then we have to make some
        %adjustments to make the matrix square before switching (see below).
        
        %If not all markers are detected (curr_num_clust<num_clust), make squared distance matrix square
        %(so that matching works below). We add additional columns of zeros so that
        %the matrix has size num_clust x num_clust.         
        if curr_num_clust<num_clust
            temp2=horzcat(temp2,zeros(num_clust,num_clust-curr_num_clust));
        end
        
        %If more clusters are detected than there are markers, we add
        %additional rows of zeros so that the matrix has size
        %curr_num_clust x curr_num_clust
        if curr_num_clust>num_clust
            temp2=vertcat(temp2,zeros(curr_num_clust-num_clust,curr_num_clust));
        end
        
        %Test all permuations of possible clusters for best assignment
        perm_mat=perms(1:max(num_clust,curr_num_clust)); %Matrix of all permutations of the number of clusters or markers (whichever is larger)
        val=zeros(1,size(perm_mat,1)); %This keeps track of the value we're trying to minimize for each possible switch (permutation)
        for p=1:size(perm_mat,1); %Loop through all permutations
            val(p)=trace(temp2(:,perm_mat(p,:))); %This is what we're trying to minimize (sum of squared distances of diagonal terms in the matrix)
        end
        
        %Get the permutation that led to the smallest sum of squared
        %distance of diagonal terms (val)
        [~,min_idx]=min(val); 
        perm_idx=perm_mat(min_idx,:);
        
        %meds (a temporary variable) lists the marker positions. If
        %some markers are missing in this frame (curr_num_clust<num_clust),
        %then add NaNs to meds for each missing marker.
        if curr_num_clust<num_clust
            meds=vertcat(meds,NaN(num_clust-curr_num_clust,3));
        end
        
        %meds2 (a temporary variable) lists the marker positions after
        %assignment.
        meds2=meds(perm_idx,:); %Switch the marker positions to get the correct assignments (according to perm_idx we found above)
        %If there are more clusters than makers, remove the final cluster
        %(which corresponds to the one that matched the previous markers
        %the least).
        if curr_num_clust>num_clust
            meds2(curr_num_clust,:)=[]; %Note to self - this only deals w/ when curr_num_clust==num_clust+1 (which is all we deal with)
        end
        
        %Update the medians matrix (which has the marker positions at all
        %times) for the current frame. This matrix has NaNs when markers are missing.
        medians(:,:,t)=meds2; 
        
        %Update the medians2 matrix (which has the marker positions at all
        %times) for the current frame. This matrix has the previous known position when markers are missing.
        medians2(:,:,t)=meds2;
        nanidx=find(isnan(medians2(:,1,t))); %Find missing markers
        medians2(nanidx,:,t)=medians2(nanidx,:,t-1); %Fill them in with value at last time point (which will be the last known position)
        
    end
    
    %Update prev_meds: The locations of the markers detected in the current frame
    %will be the previous locations of markers at the next time point
    prev_meds=meds2;
    prev_meds=prev_meds(~isnan(prev_meds(:,1)),:); %Only include the markers that were detected (not the missing ones w/ NaNs)
    
    %Update prev_num_clust. This will generally be number of clusters in
    %the current frame (which will be the number of clusters in the
    %previous frame, at the next time point). However, if more clusters
    %were detected than markers, we set prev_num_clust to be the number of
    %markers (since the extra clusters were discarded).
    prev_num_clust=min(curr_num_clust,num_clust); 
     
end

end