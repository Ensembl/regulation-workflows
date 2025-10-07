# Signal generation workflow

## Submitting Signal generation jobs

To submit the workflow, use the following command:

```
argo submit \
--namespace argo \
--serviceaccount ensreg \
--from workflowtemplate/get-signal-task-v-0.2-0
--parameter-file <JSON or YAML file with parameters>
```

Parameter file example (YAML):
```yaml
species_name: "Bos taurus"
experiment_type: "atac_seq"
skip: "0"
limit: "20"
```

## Data flow


```mermaid
graph TD
    %% Entry Point
    subgraph submitter["Signal Generation Tasks Submitter"]
        INPUT[/"Parameters:<br/>• species/assembly<br/>• experiment filters<br/>• epigenome filters<br/>• p-value-based flag<br/>• skip/limit"/]
        INPUT -->|"Filter Criteria"| A["get-signal-tasks"]
        A -->|"Query Parameters"| B[compute-get-signal-tasks-url]
        B -->|"API URL"| C[submit-get-request-from-reg-pipelines-api]
        C ==>|"Task List + Count"| D[submit-signal-tasks]
    end
    
    %% Signal Task Processing
    D -->|"Task Index"| E[submit-get-wf-task-payload-by-index]
    E ==>|"Task Payload"| F[submit-pq-values-to-bigwig]
    
    %% Work Avoidance
    F -->|"Task Payload"| G[compute-signal-task-marker-name]
    G -->|"Marker S3 Key"| H[check-if-task-marker-exists]
    F -->|"Output S3 Key"| I[check-if-bigwig-file-already-in-s3-bucket]
    H -->|"Marker Status"| J[check-work-avoidance-consistency]
    I -->|"File Status"| J
    
    %% PVC Setup
    J -->|"avoid-work == false"| K[get-signal-pvc-size]
    K -->|"PVC Size"| L[create-signal-pvc]
    L -->|"PVC Name"| M[copy-pq-value-from-s3-to-pvc]
    M -->|"PQ Values TGZ"| N[unzip-pq-values-file]
    
    %% BigWig Generation
    N -->|"PQ Values TXT"| O[genrich-pq-values-to-bigwig]
    
    subgraph bigwigGeneration["BigWig Generation from PQ-Values"]
        O -->|"p-value-based flag"| P{Use P-values<br/>or Q-values?}
        
        P ==>|"p-value-based == true"| Q["create-bed-from-p-values<br/>(awk: extract column NF-1)"]
        P ==>|"p-value-based == false"| R["create-bed-from-q-values<br/>(awk: extract column NF)"]
        
        Q -->|"Unsorted BED"| S["bedSort"]
        R -->|"Unsorted BED"| S
        
        S -->|"Sorted BED"| T["bedGraphToBigWig"]
        T ==>|"BigWig File"| U[save bigwig to S3]
    end
    
    %% Post-Processing
    U -->|"BigWig in PVC"| V[execute-compute-bigwig-file-metadata]
    V -->|"File Size + MD5"| W[execute-post-signal-file]
    
    %% Cleanup
    T -->|"BED Filename"| X[remove-intermediate-bed-file-from-pvc]
    T -->|"Sorted BED Filename"| Y[remove-intermediate-sorted-bed-file-from-pvc]
    
    U -->|"Success"| Z[update-signal-task-marker]
    U -->|"PVC Name"| AA[delete-signal-pvc]
    
    %% Styling
    classDef genomicTool fill:#4f46e5,color:#ffffff,stroke:#3730a3,stroke-width:2px
    classDef decision fill:#fef3c7,color:#92400e,stroke:#f59e0b,stroke-width:2px
    classDef workAvoidance fill:#e5e7eb,color:#111827,stroke:#6b7280,stroke-width:1px
    classDef storage fill:#fbcfe8,color:#831843,stroke:#db2777,stroke-width:1px
    classDef inputParams fill:#bae6fd,color:#0c4a6e,stroke:#0284c7,stroke-width:2px
    
    class INPUT inputParams
    class Q,R,S,T genomicTool
    class P decision
    class G,H,I,J workAvoidance
    class K,L,M,N,U,V,W,X,Y,Z,AA storage

```