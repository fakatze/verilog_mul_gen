# -*- coding: utf-8 -*-

import sys


class Expr:
    # Базовый класс, определяет представление узла в дереве выражений.


    wire = None
    done = False


# Описанные ниже классы, наследуют представления узлов у Expr
#   и получают аргументы от вложенного объекта, что позволяет пользоваться
#   данной конструкцией, без нарушения целостности структуры программы.

class ConstBit(Expr):
    # Декларация значенй 0 или 1.



    def __init__(self, v):
        assert v in (0, 1)
        self.v = v


class InBit(Expr):
    # Представление бита, входящего в дерево выражений.

    def __init__(self, xy, p):
        assert xy in ('x', 'y')
        self.xy = xy
        self.p = p


class NotBit(Expr):
    # Представление инвертора.

    def __init__(self, v):
        self.v = v


class BoothNeg(Expr):
    # Представление расчета инверсного флага в 4-х значном кодировщике Бута

    def __init__(self, pat):
        assert len(pat) == 3
        self.pat = pat


class BoothProd(Expr):
    # Представление расчета частного произведения в кодировке Бута.

    def __init__(self, pat, b):
        assert len(pat) == 3
        assert len(b) == 2
        self.pat = pat
        self.b = b


class AddBitD(Expr):
    # Представление разряда, выбираемого из сумматора.

    def __init__(self, v):
        self.v = v


class AddBitC(Expr):
    # Представление бита переноса.

    def __init__(self, v):
        self.v = v


class HalfAdd(Expr):
    # Представление полусумматора.

    def __init__(self, a, b):
        self.a = a
        self.b = b


class FullAdd(Expr):
    # Представление полного сумматора.

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c


class CarryProp(Expr):
    # Представление основного узла дерева битов переноса.

    def __init__(self, a, b):
        self.a = a
        self.b = b


class CarryMerge(Expr):
    # Представление внутреннего узла дерева битов переноса.

    def __init__(self, p0, p1):
        self.p0 = p0
        self.p1 = p1


class CarryEval(Expr):
    # Представление логики выбора бита переноса.

    def __init__(self, p, c):
        self.p = p
        self.c = c


def gen_partial_products(xvec, yvec):
    # Формируем список частичных произведений с помощью алгоритма Бута при помощи 4-Radix


    partial_products = []

    # Append zero on LSB side of xvec, sign-extend on MSB side of xvec.
    xtmp = [ConstBit(0)] + xvec + xvec[-1:]

    # Append zero on LSB side of yvec, sign-extend on MSB side of yvec.
    ytmp = [ConstBit(0)] + yvec + yvec[-1:]

    # Step through xvec, 2 bits at a time.
    for i in range(0, len(xvec), 2):

        # Select group of 3 bits from xvec (one bit overlap with last group).
        pat = xtmp[i:i + 3]

        # Add either 0, +1, +2, -1 or -2 times yvec according to Booth method.
        # Step through the bits of yvec.
        for j in range(len(yvec) + 1):

            # Use Booth encoder to choose between 0, yvec[j], yvec[j-1] or
            # inverted bits yvec[j] or yvec[j-1].
            t = BoothProd(pat, ytmp[j:j + 2])

            # Invert the MSB bit, except on first row.
            if i > 0 and j == len(yvec):
                t = NotBit(t)

            # Add result as partial product.
            partial_products.append((i + j, t))

        # For first row, sign-extend by two bits.
        # Apply sign inversion on the new MSB bit.
        if i == 0:
            partial_products.append((i + len(yvec) + 1, t))
            partial_products.append((i + len(yvec) + 2, NotBit(t)))

        # For each row except the first row, add constant 1 in the next column.
        if i > 0:
            partial_products.append((i + len(yvec) + 1, ConstBit(1)))

        # Use Booth encoder to add 1 in case of negative factor (-1 or -2).
        t = BoothNeg(pat)
        partial_products.append((i, t))

    return partial_products


def gen_dadda_tree(partial_products, nbits):
    """Generate carry save adder based on Dadda tree."""

    # Sort partial products by bit position.
    tvec = [[] for p in range(nbits)]
    for (p, b) in partial_products:
        if p < nbits:
            tvec[p].append(b)

    # Build Dadda tree.
    while any([len(t) > 3 for t in tvec]):
        # New layer.
        nvec = [[] for p in range(nbits + 1)]
        for p in range(nbits):
            t = tvec[p]
            i = 0
            while i + 2 < len(t):
                # build full adder
                a = FullAdd(t[i], t[i + 1], t[i + 2])
                nvec[p].append(AddBitD(a))
                nvec[p + 1].append(AddBitC(a))
                i += 3
            if i + 1 < len(t) and len(nvec[p]) % 3 == 2:
                # build half adder
                a = HalfAdd(t[i], t[i + 1])
                nvec[p].append(AddBitD(a))
                nvec[p + 1].append(AddBitC(a))
                i += 2
            if i < len(t):
                # pass through
                nvec[p] += t[i:]
        tvec = nvec[:nbits]

    # Last layer.
    nvec = [[] for p in range(nbits + 1)]
    for p in range(nbits):
        t = tvec[p]
        if len(t) == 3:
            # full adder
            a = FullAdd(t[0], t[1], t[2])
            nvec[p].append(AddBitD(a))
            nvec[p + 1].append(AddBitC(a))
        elif len(t) == 2 and len(nvec[p]) > 0:
            # half adder
            a = HalfAdd(t[0], t[1])
            nvec[p].append(AddBitD(a))
            nvec[p + 1].append(AddBitC(a))
        else:
            # pass through
            nvec[p] += t
    tvec = nvec[:nbits]

    # Extract remaining two rows of bits.
    avec = [(t[0] if len(t) > 0 else ConstBit(0)) for t in tvec]
    bvec = [(t[1] if len(t) > 1 else ConstBit(0)) for t in tvec]

    return (avec, bvec)


def gen_adder(avec, bvec):
    # Генерируем сумматор с ускоренным переносом.

    def carry_lookahead(pvec, cin):
        # Рекурсивно определяем направление переноса.

        if len(pvec) == 1:
            prop = pvec[0]
            cvec = [cin]
        else:
            k = (len(pvec) + 1) // 2
            (p0, c0) = carry_lookahead(pvec[:k], cin)
            ctmp = CarryEval(p0, cin)
            (p1, c1) = carry_lookahead(pvec[k:], ctmp)
            prop = CarryMerge(p0, p1)
            cvec = c0 + c1

        return (prop, cvec)

    assert len(avec) == len(bvec)

    # Определить генерацию переноса и распространения для каждого положения.
    pvec = [CarryProp(a, b) for (a, b) in zip(avec, bvec)]

    # Определить перенос для каждого положения.
    (prop, cvec) = carry_lookahead(pvec, ConstBit(0))

    # Формируем массив полных сумматоров.
    sumvec = [AddBitD(FullAdd(a, b, c))
              for (a, b, c) in zip(avec, bvec, cvec)]

    return sumvec


def gen_multiplier(xbit, ybit):
    # Формируем дерево выражений с учетом логики умножителя

    xvec = [InBit('x', p) for p in range(xbit)]
    yvec = [InBit('y', p) for p in range(ybit)]

    partial_products = gen_partial_products(xvec, yvec)
    (avec, bvec) = gen_dadda_tree(partial_products, xbit + ybit)

    zvec = gen_adder(avec, bvec)

    return zvec


def gen_netlist(node, wires, insts):
    """Generate netlist consisting of wires and component instances."""

    if node.done:
        # already processed this node
        return

    node.done = True

    if isinstance(node, ConstBit):
        # resolve during code generation
        node.wire = node
    elif isinstance(node, InBit):
        # resolve during code generation
        node.wire = node
    elif isinstance(node, NotBit):
        # create output wire
        node.wire = 'winv%d' % len(wires)
        wires.append(node.wire)
        # recurse
        gen_netlist(node.v, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, BoothNeg):
        # create output wire
        node.wire = 'wboothneg%d' % len(wires)
        wires.append(node.wire)
        # recurse
        for v in node.pat:
            gen_netlist(v, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, BoothProd):
        # create output wire
        node.wire = 'wboothprod%d' % len(wires)
        wires.append(node.wire)
        # recurse
        for v in node.pat:
            gen_netlist(v, wires, insts)
        for v in node.b:
            gen_netlist(v, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, AddBitD):
        # recurse
        gen_netlist(node.v, wires, insts)
        node.wire = node.v.wire + 'd'
    elif isinstance(node, AddBitC):
        # recurse
        gen_netlist(node.v, wires, insts)
        node.wire = node.v.wire + 'c'
    elif isinstance(node, HalfAdd):
        # create output wires
        node.wire = 'wadd%d' % len(wires)
        wires.append(node.wire + 'd')
        wires.append(node.wire + 'c')
        # recurse
        gen_netlist(node.a, wires, insts)
        gen_netlist(node.b, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, FullAdd):
        # create output wires
        node.wire = 'wadd%d' % len(wires)
        wires.append(node.wire + 'd')
        wires.append(node.wire + 'c')
        # recurse
        gen_netlist(node.a, wires, insts)
        gen_netlist(node.b, wires, insts)
        gen_netlist(node.c, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, CarryProp):
        # create output wires
        node.wire = 'wcarry%d' % len(wires)
        wires.append(node.wire + 'g')
        wires.append(node.wire + 'p')
        # recurse
        gen_netlist(node.a, wires, insts)
        gen_netlist(node.b, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, CarryMerge):
        # create output wires
        node.wire = 'wcarry%d' % len(wires)
        wires.append(node.wire + 'g')
        wires.append(node.wire + 'p')
        # recurse
        gen_netlist(node.p0, wires, insts)
        gen_netlist(node.p1, wires, insts)
        # create instance
        insts.append(node)
    elif isinstance(node, CarryEval):
        # create output wire
        node.wire = 'wcarry%d' % len(wires)
        wires.append(node.wire)
        # recurse
        gen_netlist(node.p, wires, insts)
        gen_netlist(node.c, wires, insts)
        # create instance
        insts.append(node)
    else:
        assert False


def verilog_inst(node):
    """Return (name, ports) for a given instance."""

    if isinstance(node, NotBit):
        name = 'smul_inverter'
        ports = (node.v.wire, node.wire)
    elif isinstance(node, BoothNeg):
        name = 'smul_booth_neg'
        ports = (node.pat[0].wire, node.pat[1].wire, node.pat[2].wire,
                 node.wire)
    elif isinstance(node, BoothProd):
        name = 'smul_booth_prod'
        ports = (node.pat[0].wire, node.pat[1].wire, node.pat[2].wire,
                 node.b[0].wire, node.b[1].wire,
                 node.wire)
    elif isinstance(node, HalfAdd):
        name = 'smul_half_add'
        ports = (node.a.wire, node.b.wire,
                 node.wire + 'd', node.wire + 'c')
    elif isinstance(node, FullAdd):
        name = 'smul_full_add'
        ports = (node.a.wire, node.b.wire, node.c.wire,
                 node.wire + 'd', node.wire + 'c')
    elif isinstance(node, CarryProp):
        name = 'smul_carry_prop'
        ports = (node.a.wire, node.b.wire,
                 node.wire + 'g', node.wire + 'p')
    elif isinstance(node, CarryMerge):
        name = 'smul_carry_merge'
        ports = (node.p0.wire + 'g', node.p0.wire + 'p',
                 node.p1.wire + 'g', node.p1.wire + 'p',
                 node.wire + 'g', node.wire + 'p')
    elif isinstance(node, CarryEval):
        name = 'smul_carry_eval'
        ports = (node.p.wire + 'g', node.p.wire + 'p', node.c.wire,
                 node.wire)
    else:
        assert False

    return (name, ports)


def verilog_wire(wire):
    """Resolve wire to Verilog expression string."""

    if isinstance(wire, ConstBit):
        return "1'b%d" % wire.v
    elif isinstance(wire, InBit):
        return "%sin[%d]" % (wire.xy, wire.p)
    else:
        assert isinstance(wire, str)
        return wire


def gen_verilog_lib():
    """Generate Verilog code for library components."""

    template = """

// Inverter.

module smul_inverter (
    input  wire d,
    output wire q );

assign q = ~d;

endmodule


// Half-adder.

module smul_half_add (
    input  wire x,
    input  wire y,
    output wire d,
    output wire c );

assign d = x ^ y;
assign c = x & y;

endmodule


// Full-adder.

module smul_full_add (
    input  wire x,
    input  wire y,
    input  wire z,
    output wire d,
    output wire c );

assign d = x ^ y ^ z;
assign c = (x & y) | (y & z) | (x & z);

endmodule


// Booth negative flag.

module smul_booth_neg (
    input  wire p0,
    input  wire p1,
    input  wire p2,
    output wire f );

assign f = p2 & ((~p1) | (~p0));

endmodule


// Booth partial product generator.

module smul_booth_prod (
    input  wire p0,
    input  wire p1,
    input  wire p2,
    input  wire u0,
    input  wire u1,
    output reg  y );

always @ (*)
begin
    case ({p2, p1, p0})
        3'b000  : y = 1'b0;
        3'b001  : y = u1;
        3'b010  : y = u1;
        3'b011  : y = u0;
        3'b100  : y = ~u0;
        3'b101  : y = ~u1;
        3'b110  : y = ~u1;
        default : y = 1'b0;
    endcase
end

endmodule


// Deterimine carry generate and carry propagate.

module smul_carry_prop (
    input  wire a,
    input  wire b,
    output wire g,
    output wire p );

assign g = a & b;
assign p = a ^ b;

endmodule


// Merge two carry propagation trees.

module smul_carry_merge (
    input  wire g0,
    input  wire p0,
    input  wire g1,
    input  wire p1,
    output wire g,
    output wire p );

assign g = g1 | (g0 & p1);
assign p = p0 & p1;

endmodule


// Calculate carry-out through a carry propagation tree.

module smul_carry_eval (
    input  wire g,
    input  wire p,
    input  wire cin,
    output wire cout );

assign cout = g | (p & cin);

endmodule"""

    print(template)


def main():
    print("1st bitdepth:")
    global xbit
    xbit = int(input())
    print("2nd bitdepth:")
    ybit = int(input())

    if xbit < 4 or ybit < xbit:
        print('less then 4bits!')
        sys.exit(1)
    else:
        # Начинаем формирование дерева выражений.
        zvec = gen_multiplier(xbit, ybit)

        # Генерируем связи и модули.
        wires = []
        insts = []
        for node in zvec:
            gen_netlist(node, wires, insts)
        outputs = [ node.wire for node in zvec]


        xxbit=xbit-1
        yybit=ybit-1
        zleft=ybit+xbit-1
        # Упаковываем библиотеку внутренних модулей.
        gen_verilog_lib()


        print('/*')
        # print("*",genx , "x",geny, "bit signed multiplier." )
        print("* %d x %d bit signed multiplier." % (xbit, ybit))
        print('*')
        print('*/')

        print("module smul_%d_%d (" % (xbit, ybit))
        print("    input  wire [%d:0] xin," % xxbit)
        print("    input  wire [%d:0] yin," % yybit)
        print("    output wire [%d:0] zout );" % zleft)
        print('')
        # Declare signals.
        for w in wires:
            print("wire %s;" % w)

            # Instantiate components.
        for (i, node) in enumerate(insts):
            (name, ports) = verilog_inst(node)
            print("%s u%d (" % (name, i),
            ", ".join([verilog_wire(p) for p in ports]), ");")
        print('')

            # Drive output signals.
        for (i, wire) in enumerate(outputs):
            print("assign zout[%d] = %s;" % (i, verilog_wire(wire)))

            # End module.
        print('')
        print("endmodule")


if __name__ == '__main__':
    main()
