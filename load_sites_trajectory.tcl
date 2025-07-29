# 设置参数
set total_frames 10000
set step 1000
set path "./20250727purediffusion/"   ;# 改为你的 xyz 文件目录路径（可为绝对路径）
set molid -1

# 加载第一个帧文件初始化轨迹
set filename [format "%ssites_%d.xyz" $path 0]
mol new $filename type xyz waitfor all
set molid [molinfo top]

# 加载其余帧
for {set i 1} {$i < $total_frames} {incr i} {
    set index [expr $i * $step]
    set filename [format "%ssites_%d.xyz" $path $index]
    if {[file exists $filename]} {
        puts "Loading $filename"
        mol addfile $filename type xyz waitfor all mol $molid
    } else {
        puts "WARNING: $filename not found, skipping."
    }
}

# 自动播放轨迹动画
animate goto 0
animate forward
animate play