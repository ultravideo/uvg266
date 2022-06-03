import sys
from pathlib import Path
from pprint import pprint


def get_ctx_count_and_names(cabac_source_file: Path):
    count = 0
    names = []
    with open(cabac_source_file) as file:
        for line in file:
            line = line.strip().split(";")[0]
            if line.startswith("cabac_ctx_t") and not "*" in line:
                temp = 1
                name = line.split()[1].split("[")[0]
                for k in line.split("[")[1:]:
                    temp *= int(k[:-1])
                count += temp
                if temp == 1:
                    names.append(name)
                else:
                    names.extend([f"{name}[{x}]" for x in range(temp)])
    return count, names


def main(state_file: Path, ctx_names: list, ctx_count: int = 332, ctx_size: int = 6):
    ctx_store = dict()
    e_store = set()
    was_zero_last = False
    frame_num = -1
    with open(state_file, "rb") as file:
        try:
            while True:
                type_, x, y, depth = file.read(13).decode().split()
                # Reset stored data at the beginning of the frame
                if x == '0' and y == '0' and type_ == "S":
                    if not was_zero_last:
                        frame_num += 1
                        ctx_store = dict()
                        e_store = set()
                    was_zero_last = True
                else:
                    was_zero_last = False

                ctx = file.read(ctx_count * ctx_size)
                if type_ == "S":
                    # These shouldn't happen but just to make sure everything is working as intended
                    if ctx_store.get((x, y, depth)):
                        raise RuntimeError
                    ctx_store[(x, y, depth)] = ctx
                else:
                    if (x, y, depth) in e_store:
                        raise RuntimeError
                    e_store.add((x, y, depth))
                    if (s_ctx := ctx_store[(x, y, depth)]) != ctx:
                        actual_problem = False

                        for i in range(ctx_count):
                            temp_s = s_ctx[i * ctx_size: (i + 1) * ctx_size]
                            temp_e = ctx[i * ctx_size: (i + 1) * ctx_size]
                            if temp_s != temp_e:
                                if ctx_names[i] in ignore_list:
                                    continue
                                actual_problem = True
                                print(f"MISSMATCH in {ctx_names[i]} {frame_num=} {x=} {y=} {depth=}")
                                print(
                                    f"GOT     : {int.from_bytes(temp_s[0:2], 'little')}:"
                                    f"{int.from_bytes(temp_s[2:4], 'little')} "
                                    f"rate={int.from_bytes(temp_s[4:5], 'big')}")
                                print(
                                    f"EXPECTED: {int.from_bytes(temp_e[0:2], 'little')}:"
                                    f"{int.from_bytes(temp_e[2:4], 'little')} "
                                    f"rate={int.from_bytes(temp_e[4:5], 'big')}")
                        if actual_problem:
                            exit(1)
        except ValueError:
            # EOF
            pass


if __name__ == '__main__':
    ignore_list = {"sao_type_idx_model", "sao_merge_flag_model"}
    if len(sys.argv) < 2:
        print("Usage: name of the file storing the cabac states")
        exit(1)
    path = Path(__file__) / "../.." / "src" / "cabac.h"
    print(path.resolve())
    counts, names_ = get_ctx_count_and_names(path.resolve())
    main(Path(sys.argv[1]), names_, counts)
