function [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func(t, loc, num_clust, prev_num_clust, dist_min, penalty, prev_meds, medians, medians2 )
%This function clusters the points, and assigns them to markers.

%cluster_func2 is an updated version of this function, and has detailed
%comments. This function is different because it initializes with the
%number of clusters being equal to the number of markers (rather than being
%one greater). This leads to occasional errors, but an increase in speed.
%Thus, this function is used for preliminary clustering/assignment runs
%(and later cluster_func2 is used for more precision)


if isempty(loc)
    medians(:,:,t)=NaN;
    medians2(:,:,t)=medians2(:,:,t-1);
else
    
    curr_num_clust=num_clust;
    med_dist=zeros(1,num_clust);
    while (any(med_dist<dist_min) && curr_num_clust>=1)
        
        if size(loc,1)<curr_num_clust %In case there are fewer points than clusters I'm fitting
            curr_num_clust=size(loc,1);
        end
        
        if prev_num_clust~=curr_num_clust %Potentially could use known starting positions?
            clust=kmedoids(loc,curr_num_clust);
        else
            clust=kmeans(loc,curr_num_clust,'Start',prev_meds);
        end
        
        %3. Record median of each cluster
        %And record how many are in each cluster (this might be useful for
        %determining how many are in each cluster)
        
        meds=zeros(curr_num_clust,3);
        for j=1:curr_num_clust
            if nnz(clust==j)>1
                meds(j,:)=median(loc(clust==j,:));
            else if nnz(clust==j)==1
                    meds(j,:)=loc(clust==j,:);
                else
                    meds(j,:)=NaN(1,3);
                end
            end
        end
        med_dist=pdist(meds);
        
        if any(med_dist<dist_min) || any(isnan(med_dist))
            curr_num_clust=curr_num_clust-1;
        end
        
    end
    prev_meds=meds;
    prev_meds=prev_meds(~isnan(prev_meds(:,1)),:);
    
    %4. Assign to objects (based on how close it is to previous locations - or some estimate of current locations)
    
    %     prev=medians(:,:,t-1)
    
    %If the first frame, just assign. Otherwise, match up with assignment
    %from previous frames
    if t==1
        %For now, all of the clusters have to be their in the first frame analyzed
        medians(:,:,t)=meds;
        medians2(:,:,t)=meds;
        
    else
        temp=pdist2(medians2(:,:,t-1),meds); %Distance between all current markers and markers in previous frame
        temp2=temp.^2; %Distances squared
        
        %If not all markers detected, make squared distance matrix square
        %(so that matching works below).
        if curr_num_clust<num_clust
            temp2=horzcat(temp2,zeros(num_clust,num_clust-curr_num_clust));
        end
        
        %Test all permuations of possible clusters for best assignment
        perm_mat=perms(1:num_clust); %Matrix of all permutations
        val=zeros(1,num_clust);
        for p=1:size(perm_mat,1);
            val(p)=trace(temp2(:,perm_mat(p,:)));
        end
%         if t>1290
%             t
%         end
        %Add in penalties for switching (e.g. marker 1 to marker 2)
        if curr_num_clust<num_clust
            %             temp2=horzcat(temp2,zeros(num_clust,num_clust-curr_num_clust));
            nanidx=find(isnan(medians(:,1,t-1))); %Find which marker wasn't visible last frame
            for p=1:size(perm_mat,1);
                %Find which marker isn't visible this frame
                invis_idx=find(perm_mat(p,:)==num_clust);
                if invis_idx~=nanidx
                    val(p)=val(p)+penalty;
                end
            end
        end
        
        [~,min_idx]=min(val);
        perm_idx=perm_mat(min_idx,:);
        
        if curr_num_clust<num_clust
            meds=vertcat(meds,NaN(num_clust-curr_num_clust,3));
        end
        
        
        %Include medians of visible markers for current frame
        medians(:,:,t)=meds(perm_idx,:);
        
        %Include medians of all markers for current frame
        medians2(:,:,t)=meds(perm_idx,:);
        nanidx=find(isnan(medians2(:,1,t)));
        medians2(nanidx,:,t)=medians2(nanidx,:,t-1);
        
    end
    prev_num_clust=curr_num_clust;
    
end


end

