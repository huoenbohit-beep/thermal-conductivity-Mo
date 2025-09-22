#!/bin/bash

# ���ø�·��
base="/share/home/e6a000185/now work/Mo_eam.fs_tiaocan/bin"
src_dir="$base/pe"
out_dir="$base"  # ���Ŀ¼���� pe ͬ��

# ���������¶�Ŀ¼��300K �� 2800K��
for temp_path in "$src_dir"/*K; do
    T=$(basename "$temp_path")  # ��ȡ�¶��ַ��������� 400K

    # �������� pe_Mo_${T}_cdpa* ��Ŀ¼
    for folder in "$temp_path"/pe_Mo_${T}_cdpa*; do
        if [[ -d "$folder" ]]; then
            # ��ȡ cdpa ��ֵ
            cdpa=$(basename "$folder" | sed 's/.*cdpa//')

            # ���� .lammpstrj �ļ�·��
            lammpstrj_file="$folder/Mo_${T}_cdpa${cdpa}.lammpstrj"
            output_file="$out_dir/Mo-${T}-cdpa${cdpa}.pe.dat"

            # ����ļ�����
            if [[ -f "$lammpstrj_file" ]]; then
                echo "Extracting: $lammpstrj_file -> $output_file"
                awk '
                    /ITEM: NUMBER OF ATOMS/ {getline; print $1 > "'"$output_file"'"}
                    /ITEM: ATOMS/ {
                        while (getline && NF > 1)
                            print $NF >> "'"$output_file"'"
                    }
                ' "$lammpstrj_file"
            else
                echo "Warning: $lammpstrj_file not found, skipping..."
            fi
        fi
    done
done

echo "Done. All potential energy data extracted to $out_dir"
