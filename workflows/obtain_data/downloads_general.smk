import os

# variables
ROOT = os.path.dirname(os.path.dirname(os.getcwd()))
RAW_DIR = os.path.join(ROOT,"data","raw")
SUPPORT_DIR = os.path.join(ROOT,"support")
SAVE_PARAMS = {"sep":"\t", "index":False, "compression":"gzip"}


##### RULES #####
rule all:
    input:
        # download SplicingLore
        os.path.join(RAW_DIR,"SplicingLore")
        
        
rule download_splicinglore:
    input:
        urls = os.path.join(SUPPORT_DIR, "SplicingLore-urls.txt")
    output:
        folder = directory(os.path.join(RAW_DIR,"SplicingLore"))
    shell:
        """
        set -eo pipefail
        
        mkdir -p {output.folder}
        
        while IFS= read -r url; do
        
            # Download the URL using wget
            wget \
                --user-agent="Chrome" \
                --no-clobber \
                --no-check-certificate \
                $url \
                --directory-prefix {output.folder}
                
        done < {input.urls}
        
        echo "Done!"
        """
