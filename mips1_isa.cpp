/**
 * @file      mips1_isa.cpp
 * @author    Sandro Rigo
 *            Marcus Bartholomeu
 *            Alexandro Baldassin (acasm information)
 *
 *            The ArchC Team
 *            http://www.archc.org/
 *
 *            Computer Systems Laboratory (LSC)
 *            IC-UNICAMP
 *            http://www.lsc.ic.unicamp.br/
 *
 * @version   1.0
 * @date      Mon, 19 Jun 2006 15:50:52 -0300
 * 
 * @brief     The ArchC i8051 functional model.
 * 
 * @attention Copyright (C) 2002-2006 --- The ArchC Team
 * 
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License 
 * along with this program; if not, write to the Free Software 
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include  "mips1_isa.H"
#include  "mips1_isa_init.cpp"
#include  "mips1_bhv_macros.H"
#include <ctype.h>


//If you want debug information for this model, uncomment next line
//#define DEBUG_MODEL
#include "ac_debug_model.H"

int TestaDataHazard();
int TestaControlHazard();
void Empty();

typedef enum{
	iempty, ilb, ilbu, ilh, ilhu, ilw, ilwl, ilwr, isb, ish, isw, iswl, iswr, iaddi, iaddiu, islti, isltiu, iandi, iori, ixori, ilui, iadd, iaddu, isub, isubu, islt, isltu, iinstr_and, iinstr_or, iinstr_xor, iinstr_nor, inop, isll, isrl, isra, isllv, isrlv, israv, imult, imultu, idiv, idivu, imfhi, imthi, imflo, imtlo, ij, ijal, ijr, ijalr, ibeq, ibne, iblez, ibgtz, ibltz, ibgez, ibltzal, ibgezal, isys_call, iinstr_break
}TipoInstrucao;

typedef struct{
	TipoInstrucao instrucao;
	int rs;
	int rt;
	int rd;
	char type;
}Instrucao;

//!User defined macros to reference registers.
#define Ra 31
#define Sp 29



/* Defines sobre o modelo */
/* estagios do pipeline */
int N_STAGES;
/* bits do branch predictor (1 = 1 bit, 0 = 2 bits) */
int BR_PR;
/* Superscalar? 1 = Sim, 0 = Nao */
int SUPERSCALAR;

/* Variaveis globais*/
/* contadores */
int ciclos;
int bolhas;

/* 'bit' do branch predictor: */
/* se for 1bit: 0: not taken, 1: taken */
/* se for 2bit: 0 e 1: not taken, 2 e 3: taken */
int historicoPredictor; 
/* referencia: 1 se branch foi realizado. 0 caso contrario */
int taken;

/*tipo de instrucao - R, I ou J */
char InsType;

/* Snapshot do Pipeline */
Instrucao* hist;

/* Acrescenta nova instrucao,e respectivos registradores, no vetor hist[] */
/* Tambem verifica hazards e outras otimizacoes, atualizando os contadores */
void Put(Instrucao novaInstrucao){
	int i = 0,aux2;
	for (i = 0; i < N_STAGES-1; i++){
		hist[i + 1] = hist[i];
	}
	hist[0] = novaInstrucao;
	hist[0].type = InsType;
	/* TestaControlHazard() apenas se for jump ou branch - afeta instrucoes futuras */
	if(hist[0].instrucao >= ij && hist[0].instrucao <= ibgezal) {
		ciclos += TestaControlHazard();
	}
	/* verifica se houve economia de ciclo com SuperScalar - afeta instrucao passada */
	/* verifica tambem os data hazards - depende se for superscalar ou nao */
	ciclos += TestaDataHazard();
	/* soma ciclos das bolhas geradas, e 1 ciclo da propria instrucao */
	ciclos ++;
}

/* retorna true se houver data hazard */
int dependencia(Instrucao t1, Instrucao t2) {
	if(t2.type == 'X') return false;
	if(t1.type == 'R') {
		if(t2.type == 'R') 
			return (t2.rd == t1.rs || t2.rd == t1.rt);
		else if(t2.type == 'I')
			return (t2.rt == t1.rs || t2.rt == t1.rt);
	} else if(t1.type == 'I' || t1.type == 'J') {
		if(t2.type == 'I')
			return (t1.rs == t2.rt);
		else if(t2.type == 'R')
			return (t1.rs == t2.rd);
	}
	printf("Erro na checagem de dependencia!\n");
	return false;
}

/* insere n bolhas na frente da instrucao na primeira posicao do vetor historico */
/* para esvaziar o pipeline (n == N_STAGES), usar Empty() */
void addbolhas(int n) {
	int i;
	
	for(i=N_STAGES-1;i>n;i--) {
		hist[i] = hist[i-n];
	}
	for(i=1;i<n;i++) {
		hist[i].type = 'X';
		hist[i].rd = -1;
		hist[i].rs = -1;
		hist[i].rt = -1;
		hist[i].instrucao = iempty;
	}
	bolhas += n;
	return;
}

int TestaDataHazard() {
	if(false) {
		switch(N_STAGES) {
			case 5: {
				
			}
			case 7: {
			
			}
			case 13: {
			
			}
		}
	}
	else {
		switch(N_STAGES) {
			case 5: {
				/*verifica se rd das duas instrucoes anteriores sao	rt ou rs da atual */
				if(dependencia(hist[0],hist[1])) { addbolhas(2); return 2; }
				else if(dependencia(hist[0],hist[2])) { addbolhas(1); return 1; }
				break;
			}
			case 7: {
			    /*Instruction Fetch tem um estagio a mais, e MEM tambem tem. Eh preciso 
			      checar 3 instrucoes anteriores */
				if(dependencia(hist[0],hist[1])) { addbolhas(3); return 3; }
				else if(dependencia(hist[0],hist[2])) { addbolhas(2); return 2; }
				else if(dependencia(hist[0],hist[3])) { addbolhas(1); return 1; }
				break;
			}
			case 13: {
				/*Instruction Fetch:2 estagios, ID: 5 estagios, REG acessado no estagio 8
				 e escrito no estagio 13 */
				 if(dependencia(hist[0],hist[1])) { addbolhas(5); return 5; }
				 if(dependencia(hist[0],hist[2])) { addbolhas(4); return 4; }
				 if(dependencia(hist[0],hist[3])) { addbolhas(3); return 3; }
				 if(dependencia(hist[0],hist[4])) { addbolhas(2); return 2; }
				 if(dependencia(hist[0],hist[5])) { addbolhas(1); return 1; }
				 break;
			}
		}
	}
	return 0;
}

/* Verifica o status do branch predictor e se houve ou não o branch, e computa
ciclos a mais (bolha) ou nao. Se nao houve bolha, retorna 0 */
int TestaControlHazard(){

	/* se for algum Jump, bolhas geradas sao os estagios anteriores a MEM, que eh
	   quando o PC eh atualizado */
	if (hist[0].instrucao >= ij && hist[0].instrucao <= ijalr) {
		Empty();
		bolhas += N_STAGES-2;
		return N_STAGES-2;
	}
	
	/*se for algum branch */
	if(hist[0].instrucao >= ibeq && hist[0].instrucao <= ibgezal) {
		/*se predictor for de 1 bit */
		if (BR_PR) {
			/* se predictor acertou, nao ha bolhas extras*/
			if(taken == historicoPredictor) return 0;
			/*se errou, tem que inserir bolhas e atualizar bit historico*/
			else {
				Empty();
				historicoPredictor = taken;
				/* N_STAGES (encher pipeline novamente) */
				bolhas += N_STAGES;
				return N_STAGES;
			}
		}
		/* se predictor for de 2 bits */
		else {
			if(taken && (historicoPredictor == 2 || historicoPredictor == 3)) {
				historicoPredictor = 3;
				return 0;
			}
			else if(taken && historicoPredictor) {
				historicoPredictor = 3;
				bolhas += N_STAGES;
				Empty();
				return N_STAGES;
			}
			else if(taken && !historicoPredictor) {
				historicoPredictor = 1;
				bolhas += N_STAGES;
				Empty();
				return N_STAGES;
			}
			else if(!taken && (historicoPredictor == 0 || historicoPredictor == 1)) {
				historicoPredictor = 0;
				return 0;
			}
			else if(!taken && historicoPredictor == 3) {
				historicoPredictor = 2;
				bolhas += N_STAGES;
				Empty();
				return N_STAGES;
			}
			else if(!taken && historicoPredictor == 2) {
				historicoPredictor = 0;
				bolhas += N_STAGES;
				Empty();
				return N_STAGES;
			}
		}
		printf("Erro no Branch Predictor!\n");
	}
	return 0;
}

/*Funcao que esvazia vetor historico
*/
void Empty()
{
	int i;
	for (i = 0; i < N_STAGES; i++)
	{
		hist[i].instrucao = iempty;
		hist[i].rs = -1;
		hist[i].rt = -1;
		hist[i].rd = -1;
		hist[i].type = 'X';
	}
	return;
}


// 'using namespace' statement to allow access to all
// mips1-specific datatypes
using namespace mips1_parms;

//!Generic instruction behavior method.
void ac_behavior( instruction )
{ 
  dbg_printf("----- PC=%#x ----- %lld\n", (int) ac_pc, ac_instr_counter);
  //  dbg_printf("----- PC=%#x NPC=%#x ----- %lld\n", (int) ac_pc, (int)npc, ac_instr_counter);
#ifndef NO_NEED_PC_UPDATE
  ac_pc = npc;
  npc = ac_pc + 4;
#endif
};
 
//! Instruction Format behavior methods.
void ac_behavior( Type_R ){ InsType = 'R'; }
void ac_behavior( Type_I ){ InsType = 'I'; }
void ac_behavior( Type_J ){ InsType = 'J'; }
 
//!Behavior called before starting simulation
void ac_behavior(begin)
{
  dbg_printf("@@@ begin behavior @@@\n");
  RB[0] = 0;
  npc = ac_pc + 4;

  // Is is not required by the architecture, but makes debug really easier
  for (int regNum = 0; regNum < 32; regNum ++)
    RB[regNum] = 0;
  hi = 0;
  lo = 0;
  
  int aux;
  char c;
  printf("Numero de estagios do pipeline: ");
  scanf("%d",&N_STAGES);
  printf("Bits do Branch Prediction: ");
  scanf("%d",&aux);
  getchar();
  if(aux==2) BR_PR = 0;
  else BR_PR = 1;
  printf("Superscalar? (s/n) :");
  scanf(" %c",&c);
  if(toupper(c) == 'S') SUPERSCALAR = 1;
  else SUPERSCALAR = 0;
  
  hist = (Instrucao *)malloc(N_STAGES * sizeof(Instrucao));
  Empty();
  ciclos = N_STAGES-1;
  bolhas = 0;
  historicoPredictor = 0;
}

//!Behavior called after finishing simulation
void ac_behavior(end)
{
  int bits;
  
  dbg_printf("@@@ end behavior @@@\n");
  printf("Ciclos: %d\nBolhas: %d\n",ciclos,bolhas);
  printf("Configuracao:\n   Estagios de Pipeline: %d\n",N_STAGES);
  if(BR_PR) bits=1;
  else bits=2;
  printf("   Branch Predictor: %d bits\n",bits);
  if(SUPERSCALAR) printf("   Superscalar Architecture\n");
  else printf("   Scalar Architecture\n");
}


//!Instruction lb behavior method.
void ac_behavior( lb )
{
  char byte;
  
  Instrucao t;
  t.instrucao = ilb;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("lb r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = DM.read_byte(RB[rs]+ imm);
  RB[rt] = (ac_Sword)byte ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lbu behavior method.
void ac_behavior( lbu )
{
  unsigned char byte;
  
  Instrucao t;
  t.instrucao = ilbu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("lbu r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = DM.read_byte(RB[rs]+ imm);
  RB[rt] = byte ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lh behavior method.
void ac_behavior( lh )
{
  short int half;
  
  Instrucao t;
  t.instrucao = ilh;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("lh r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  half = DM.read_half(RB[rs]+ imm);
  RB[rt] = (ac_Sword)half ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lhu behavior method.
void ac_behavior( lhu )
{
  unsigned short int  half;
  
  Instrucao t;
  t.instrucao = ilhu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  half = DM.read_half(RB[rs]+ imm);
  RB[rt] = half ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lw behavior method.
void ac_behavior( lw )
{
  Instrucao t;
  t.instrucao = ilw;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("lw r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  RB[rt] = DM.read(RB[rs]+ imm);
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lwl behavior method.
void ac_behavior( lwl )
{
  dbg_printf("lwl r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;
  
  Instrucao t;
  t.instrucao = ilwl;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);

  addr = RB[rs] + imm;
  offset = (addr & 0x3) * 8;
  data = DM.read(addr & 0xFFFFFFFC);
  data <<= offset;
  data |= RB[rt] & ((1<<offset)-1);
  RB[rt] = data;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lwr behavior method.
void ac_behavior( lwr )
{
  dbg_printf("lwr r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;
  
  Instrucao t;
  t.instrucao = ilwr;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);

  addr = RB[rs] + imm;
  offset = (3 - (addr & 0x3)) * 8;
  data = DM.read(addr & 0xFFFFFFFC);
  data >>= offset;
  data |= RB[rt] & (0xFFFFFFFF << (32-offset));
  RB[rt] = data;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction sb behavior method.
void ac_behavior( sb )
{
  unsigned char byte;
  
  Instrucao t;
  t.instrucao = isb;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("sb r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  byte = RB[rt] & 0xFF;
  DM.write_byte(RB[rs] + imm, byte);
  dbg_printf("Result = %#x\n", (int) byte);
};

//!Instruction sh behavior method.
void ac_behavior( sh )
{
  unsigned short int half;
  
  Instrucao t;
  t.instrucao = ish;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("sh r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  half = RB[rt] & 0xFFFF;
  DM.write_half(RB[rs] + imm, half);
  dbg_printf("Result = %#x\n", (int) half);
};

//!Instruction sw behavior method.
void ac_behavior( sw )
{
  Instrucao t;
  t.instrucao = isw;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("sw r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  DM.write(RB[rs] + imm, RB[rt]);
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction swl behavior method.
void ac_behavior( swl )
{
  dbg_printf("swl r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;
  
  Instrucao t;
  t.instrucao = iswl;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);

  addr = RB[rs] + imm;
  offset = (addr & 0x3) * 8;
  data = RB[rt];
  data >>= offset;
  data |= DM.read(addr & 0xFFFFFFFC) & (0xFFFFFFFF << (32-offset));
  DM.write(addr & 0xFFFFFFFC, data);
  dbg_printf("Result = %#x\n", data);
};

//!Instruction swr behavior method.
void ac_behavior( swr )
{
  dbg_printf("swr r%d, %d(r%d)\n", rt, imm & 0xFFFF, rs);
  unsigned int addr, offset;
  ac_Uword data;
  
  Instrucao t;
  t.instrucao = iswr;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);

  addr = RB[rs] + imm;
  offset = (3 - (addr & 0x3)) * 8;
  data = RB[rt];
  data <<= offset;
  data |= DM.read(addr & 0xFFFFFFFC) & ((1<<offset)-1);
  DM.write(addr & 0xFFFFFFFC, data);
  dbg_printf("Result = %#x\n", data);
};

//!Instruction addi behavior method.
void ac_behavior( addi )
{
  Instrucao t;
  t.instrucao = iaddi;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("addi r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] + imm;
  dbg_printf("Result = %#x\n", RB[rt]);
  //Test overflow
  if ( ((RB[rs] & 0x80000000) == (imm & 0x80000000)) &&
       ((imm & 0x80000000) != (RB[rt] & 0x80000000)) ) {
    fprintf(stderr, "EXCEPTION(addi): integer overflow.\n"); exit(EXIT_FAILURE);
  }
};

//!Instruction addiu behavior method.
void ac_behavior( addiu )
{
  Instrucao t;
  t.instrucao = iaddiu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("addiu r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] + imm;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction slti behavior method.
void ac_behavior( slti )
{
  Instrucao t;
  t.instrucao = islti;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("slti r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Set the RD if RS< IMM
  if( (ac_Sword) RB[rs] < (ac_Sword) imm )
    RB[rt] = 1;
  // Else reset RD
  else
    RB[rt] = 0;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction sltiu behavior method.
void ac_behavior( sltiu )
{
  Instrucao t;
  t.instrucao = isltiu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("sltiu r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Set the RD if RS< IMM
  if( (ac_Uword) RB[rs] < (ac_Uword) imm )
    RB[rt] = 1;
  // Else reset RD
  else
    RB[rt] = 0;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction andi behavior method.
void ac_behavior( andi )
{	
  Instrucao t;
  t.instrucao = iandi;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("andi r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] & (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction ori behavior method.
void ac_behavior( ori )
{	
  Instrucao t;
  t.instrucao = iori;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("ori r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] | (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction xori behavior method.
void ac_behavior( xori )
{	
  Instrucao t;
  t.instrucao = ixori;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("xori r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  RB[rt] = RB[rs] ^ (imm & 0xFFFF) ;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction lui behavior method.
void ac_behavior( lui )
{	
  Instrucao t;
  t.instrucao = ilui;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("lui r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  // Load a constant in the upper 16 bits of a register
  // To achieve the desired behaviour, the constant was shifted 16 bits left
  // and moved to the target register ( rt )
  RB[rt] = imm << 16;
  dbg_printf("Result = %#x\n", RB[rt]);
};

//!Instruction add behavior method.
void ac_behavior( add )
{
  Instrucao t;
  t.instrucao = iadd;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("add r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] + RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
  //Test overflow
  if ( ((RB[rs] & 0x80000000) == (RB[rd] & 0x80000000)) &&
       ((RB[rd] & 0x80000000) != (RB[rt] & 0x80000000)) ) {
    fprintf(stderr, "EXCEPTION(add): integer overflow.\n"); exit(EXIT_FAILURE);
  }
};

//!Instruction addu behavior method.
void ac_behavior( addu )
{
  Instrucao t;
  t.instrucao = iaddu;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("addu r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] + RB[rt];
  //cout << "  RS: " << (unsigned int)RB[rs] << " RT: " << (unsigned int)RB[rt] << endl;
  //cout << "  Result =  " <<  (unsigned int)RB[rd] <<endl;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sub behavior method.
void ac_behavior( sub )
{
  Instrucao t;
  t.instrucao = isub;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("sub r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] - RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
  //TODO: test integer overflow exception for sub
};

//!Instruction subu behavior method.
void ac_behavior( subu )
{
  Instrucao t;
  t.instrucao = isubu;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("subu r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] - RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction slt behavior method.
void ac_behavior( slt )
{	
  Instrucao t;
  t.instrucao = islt;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("slt r%d, r%d, r%d\n", rd, rs, rt);
  // Set the RD if RS< RT
  if( (ac_Sword) RB[rs] < (ac_Sword) RB[rt] )
    RB[rd] = 1;
  // Else reset RD
  else
    RB[rd] = 0;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sltu behavior method.
void ac_behavior( sltu )
{
  Instrucao t;
  t.instrucao = isltu;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("sltu r%d, r%d, r%d\n", rd, rs, rt);
  // Set the RD if RS < RT
  if( RB[rs] < RB[rt] )
    RB[rd] = 1;
  // Else reset RD
  else
    RB[rd] = 0;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_and behavior method.
void ac_behavior( instr_and )
{
  Instrucao t;
  t.instrucao = iinstr_and;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("instr_and r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] & RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_or behavior method.
void ac_behavior( instr_or )
{
  Instrucao t;
  t.instrucao = iinstr_or;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("instr_or r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] | RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_xor behavior method.
void ac_behavior( instr_xor )
{
  Instrucao t;
  t.instrucao = iinstr_xor;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("instr_xor r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = RB[rs] ^ RB[rt];
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction instr_nor behavior method.
void ac_behavior( instr_nor )
{
  Instrucao t;
  t.instrucao = iinstr_nor;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("nor r%d, r%d, r%d\n", rd, rs, rt);
  RB[rd] = ~(RB[rs] | RB[rt]);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction nop behavior method.
void ac_behavior( nop )
{  
  Instrucao t;
  t.instrucao = inop;
  t.rt = -1;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  dbg_printf("nop\n");
};

//!Instruction sll behavior method.
void ac_behavior( sll )
{  
  Instrucao t;
  t.instrucao = isll;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("sll r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = RB[rt] << shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srl behavior method.
void ac_behavior( srl )
{
  Instrucao t;
  t.instrucao = isrl;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("srl r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = RB[rt] >> shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sra behavior method.
void ac_behavior( sra )
{
  Instrucao t;
  t.instrucao = isra;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("sra r%d, r%d, %d\n", rd, rs, shamt);
  RB[rd] = (ac_Sword) RB[rt] >> shamt;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction sllv behavior method.
void ac_behavior( sllv )
{
  Instrucao t;
  t.instrucao = isllv;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("sllv r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = RB[rt] << (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srlv behavior method.
void ac_behavior( srlv )
{
  Instrucao t;
  t.instrucao = isrlv;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("srlv r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = RB[rt] >> (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction srav behavior method.
void ac_behavior( srav )
{
  Instrucao t;
  t.instrucao = israv;
  t.rt = rt;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("srav r%d, r%d, r%d\n", rd, rt, rs);
  RB[rd] = (ac_Sword) RB[rt] >> (RB[rs] & 0x1F);
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mult behavior method.
void ac_behavior( mult )
{
  Instrucao t;
  t.instrucao = imult;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("mult r%d, r%d\n", rs, rt);

  long long result;
  int half_result;

  result = (ac_Sword) RB[rs];
  result *= (ac_Sword) RB[rt];

  half_result = (result & 0xFFFFFFFF);
  // Register LO receives 32 less significant bits
  lo = half_result;

  half_result = ((result >> 32) & 0xFFFFFFFF);
  // Register HI receives 32 most significant bits
  hi = half_result ;

  dbg_printf("Result = %#llx\n", result);
};

//!Instruction multu behavior method.
void ac_behavior( multu )
{
  Instrucao t;
  t.instrucao = imultu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("multu r%d, r%d\n", rs, rt);

  unsigned long long result;
  unsigned int half_result;

  result  = RB[rs];
  result *= RB[rt];

  half_result = (result & 0xFFFFFFFF);
  // Register LO receives 32 less significant bits
  lo = half_result;

  half_result = ((result>>32) & 0xFFFFFFFF);
  // Register HI receives 32 most significant bits
  hi = half_result ;

  dbg_printf("Result = %#llx\n", result);
};

//!Instruction div behavior method.
void ac_behavior( div )
{
  Instrucao t;
  t.instrucao = idiv;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("div r%d, r%d\n", rs, rt);
  // Register LO receives quotient
  lo = (ac_Sword) RB[rs] / (ac_Sword) RB[rt];
  // Register HI receives remainder
  hi = (ac_Sword) RB[rs] % (ac_Sword) RB[rt];
};

//!Instruction divu behavior method.
void ac_behavior( divu )
{
  Instrucao t;
  t.instrucao = idivu;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("divu r%d, r%d\n", rs, rt);
  // Register LO receives quotient
  lo = RB[rs] / RB[rt];
  // Register HI receives remainder
  hi = RB[rs] % RB[rt];
};

//!Instruction mfhi behavior method.
void ac_behavior( mfhi )
{
  Instrucao t;
  t.instrucao = imfhi;
  t.rt = -1;
  t.rs = -1;
  t.rd = rd;
  Put(t);
  
  dbg_printf("mfhi r%d\n", rd);
  RB[rd] = hi;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mthi behavior method.
void ac_behavior( mthi )
{
  Instrucao t;
  t.instrucao = imthi;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("mthi r%d\n", rs);
  hi = RB[rs];
  dbg_printf("Result = %#x\n", hi);
};

//!Instruction mflo behavior method.
void ac_behavior( mflo )
{
  Instrucao t;
  t.instrucao = imflo;
  t.rt = -1;
  t.rs = -1;
  t.rd = rd;
  Put(t);
  
  dbg_printf("mflo r%d\n", rd);
  RB[rd] = lo;
  dbg_printf("Result = %#x\n", RB[rd]);
};

//!Instruction mtlo behavior method.
void ac_behavior( mtlo )
{
  Instrucao t;
  t.instrucao = imtlo;
  t.rt = rs;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  dbg_printf("mtlo r%d\n", rs);
  lo = RB[rs];
  dbg_printf("Result = %#x\n", lo);
};

//!Instruction j behavior method.
void ac_behavior( j )
{
  Instrucao t;
  t.instrucao = ij;
  t.rt = -1;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  dbg_printf("j %d\n", addr);
  addr = addr << 2;
#ifndef NO_NEED_PC_UPDATE
  npc =  (ac_pc & 0xF0000000) | addr;
#endif 
  dbg_printf("Target = %#x\n", (ac_pc & 0xF0000000) | addr );
};

//!Instruction jal behavior method.
void ac_behavior( jal )
{
  Instrucao t;
  t.instrucao = ijal;
  t.rt = -1;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  dbg_printf("jal %d\n", addr);
  // Save the value of PC + 8 (return address) in $ra ($31) and
  // jump to the address given by PC(31...28)||(addr<<2)
  // It must also flush the instructions that were loaded into the pipeline
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
	
  addr = addr << 2;
#ifndef NO_NEED_PC_UPDATE
  npc = (ac_pc & 0xF0000000) | addr;
#endif 
	
  dbg_printf("Target = %#x\n", (ac_pc & 0xF0000000) | addr );
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction jr behavior method.
void ac_behavior( jr )
{
  Instrucao t;
  t.instrucao = ijr;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
  dbg_printf("jr r%d\n", rs);
  // Jump to the address stored on the register reg[RS]
  // It must also flush the instructions that were loaded into the pipeline
#ifndef NO_NEED_PC_UPDATE
  npc = RB[rs], 1;
#endif 
  dbg_printf("Target = %#x\n", RB[rs]);
};

//!Instruction jalr behavior method.
void ac_behavior( jalr )
{
  Instrucao t;
  t.instrucao = ijalr;
  t.rt = -1;
  t.rs = rs;
  t.rd = rd;
  Put(t);
  
  dbg_printf("jalr r%d, r%d\n", rd, rs);
  // Save the value of PC + 8(return address) in rd and
  // jump to the address given by [rs]

#ifndef NO_NEED_PC_UPDATE
  npc = RB[rs], 1;
#endif 
  dbg_printf("Target = %#x\n", RB[rs]);

  if( rd == 0 )  //If rd is not defined use default
    rd = Ra;
  RB[rd] = ac_pc+4;
  dbg_printf("Return = %#x\n", ac_pc+4);
};

//!Instruction beq behavior method.
void ac_behavior( beq )
{
  
  
  dbg_printf("beq r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  if( RB[rs] == RB[rt] ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }
  else taken = 0;
  
  Instrucao t;
  t.instrucao = ibeq;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction bne behavior method.
void ac_behavior( bne )
{	
  
  
  dbg_printf("bne r%d, r%d, %d\n", rt, rs, imm & 0xFFFF);
  if( RB[rs] != RB[rt] ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	
  else taken = 0;
  
  Instrucao t;
  t.instrucao = ibne;
  t.rt = rt;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction blez behavior method.
void ac_behavior( blez )
{
  dbg_printf("blez r%d, %d\n", rs, imm & 0xFFFF);
  if( (RB[rs] == 0 ) || (RB[rs]&0x80000000 ) ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2), 1;
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	
  else taken = 0;
  
  Instrucao t;
  t.instrucao = iblez;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction bgtz behavior method.
void ac_behavior( bgtz )
{ 
  dbg_printf("bgtz r%d, %d\n", rs, imm & 0xFFFF);
  if( !(RB[rs] & 0x80000000) && (RB[rs]!=0) ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	
  else taken = 0;
  
  Instrucao t;
  t.instrucao = ibgtz;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction bltz behavior method.
void ac_behavior( bltz )
{ 
  dbg_printf("bltz r%d, %d\n", rs, imm & 0xFFFF);
  if( RB[rs] & 0x80000000 ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }
  else taken = 0;
  
  Instrucao t;
  t.instrucao = ibltz;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
  
};

//!Instruction bgez behavior method.
void ac_behavior( bgez )
{
  
  
  dbg_printf("bgez r%d, %d\n", rs, imm & 0xFFFF);
  if( !(RB[rs] & 0x80000000) ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	else taken = 0;
  
  Instrucao t;
  t.instrucao = ibgez;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction bltzal behavior method.
void ac_behavior( bltzal )
{
  
  
  dbg_printf("bltzal r%d, %d\n", rs, imm & 0xFFFF);
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
  if( RB[rs] & 0x80000000 ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	else taken = 0;
  dbg_printf("Return = %#x\n", ac_pc+4);
  
  Instrucao t;
  t.instrucao = ibltzal;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction bgezal behavior method.
void ac_behavior( bgezal )
{
  
  
  dbg_printf("bgezal r%d, %d\n", rs, imm & 0xFFFF);
  RB[Ra] = ac_pc+4; //ac_pc is pc+4, we need pc+8
  if( !(RB[rs] & 0x80000000) ){
    taken = 1;
#ifndef NO_NEED_PC_UPDATE
    npc = ac_pc + (imm<<2);
#endif 
    dbg_printf("Taken to %#x\n", ac_pc + (imm<<2));
  }	else taken = 0;
  dbg_printf("Return = %#x\n", ac_pc+4);
  
  Instrucao t;
  t.instrucao = ibgezal;
  t.rt = -1;
  t.rs = rs;
  t.rd = -1;
  Put(t);
};

//!Instruction sys_call behavior method.
void ac_behavior( sys_call )
{
  Instrucao t;
  t.instrucao = isys_call;
  t.rt = -1;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  dbg_printf("syscall\n");
  stop();
}

//!Instruction instr_break behavior method.
void ac_behavior( instr_break )
{
  Instrucao t;
  t.instrucao = iinstr_break;
  t.rt = -1;
  t.rs = -1;
  t.rd = -1;
  Put(t);
  
  fprintf(stderr, "instr_break behavior not implemented.\n"); 
  exit(EXIT_FAILURE);
}
