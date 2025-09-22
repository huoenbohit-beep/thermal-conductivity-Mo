#!/bin/bash

# 设置根路径
base="/share/home/e6a000185/now work/Mo_eam.fs_tiaocan/bin"
src_dir="$base/pe"
out_dir="$base"  # 输出目录，和 pe 同级

# 遍历所有温度目录（300K 到 2800K）
for temp_path in "$src_dir"/*K; do
    T=$(basename "$temp_path")  # 提取温度字符串，比如 400K

    # 遍历所有 pe_Mo_${T}_cdpa* 子目录
    for folder in "$temp_path"/pe_Mo_${T}_cdpa*; do
        if [[ -d "$folder" ]]; then
            # 提取 cdpa 数值
            cdpa=$(basename "$folder" | sed 's/.*cdpa//')

            # 构造 .lammpstrj 文件路径
            lammpstrj_file="$folder/Mo_${T}_cdpa${cdpa}.lammpstrj"
            output_file="$out_dir/Mo-${T}-cdpa${cdpa}.pe.dat"

            # 检查文件存在
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
