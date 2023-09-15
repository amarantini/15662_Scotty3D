#include "test.h"
#include "rasterizer/pipeline.h"
#include "rasterizer/programs.h"

#include <limits>
#include <iomanip>
#include <algorithm>
#include <unordered_set>

using TestPipeline = Pipeline< PrimitiveType::Lines, Programs::Lambertian, Pipeline_Blend_Replace | Pipeline_Depth_Less | Pipeline_Interp_Flat >;

namespace std {
	template< >
	struct hash< Vec2 > {
		size_t operator()(const Vec2 &v) const {
			static hash< float > hf;
			size_t x = hf(v.x);
			size_t y = hf(v.y);
			return x ^ (y << 16) ^ (y >> (sizeof(y)*8-16));
		}
	};
}

namespace extra_test {
	//check that line produces exactly the listed fragments:
	void check_line_covers(std::string const &desc, std::vector< Vec2 > const &line_strip, std::unordered_set< Vec2 > const &expected) {

		std::unordered_set< Vec2 > got;
		for (uint32_t i = 0; i + 1 < line_strip.size(); ++i) {
			TestPipeline::ClippedVertex a,b;
			a.fb_position = Vec3(line_strip[i].x, line_strip[i].y, 0.25f);
			a.inv_w = 1.0f;
			a.attributes.fill(1.0f);
			b.fb_position = Vec3(line_strip[i+1].x, line_strip[i+1].y, 0.75f);
			b.inv_w = 1.0f;
			b.attributes.fill(2.0f);
			TestPipeline::rasterize_line(a, b, [&](TestPipeline::Fragment const &frag){
				got.emplace(frag.fb_position.x, frag.fb_position.y);
			});
		}

		std::vector< std::string > raster;

		raster.emplace_back(".");

		uint32_t out_of_raster = 0;

		auto draw = [&raster,&out_of_raster](Vec2 const &px, char c) {
			int32_t x = int32_t(std::floor(px.x));
			int32_t y = int32_t(std::floor(px.y));

			if (x < 0 || y < 0 || x > 10 || y > 10) {
				++out_of_raster;
				return;
			}

			if (uint32_t(y) >= raster.size()) {
				raster.resize(y+1, "");
			}
			if (uint32_t(x) >= raster[y].size()) {
				raster[y].resize(x+1, '.');
			}
			raster[y][x] = c;
		};

		uint32_t matched = 0;
		uint32_t missed = 0;
		uint32_t extra = 0;

		for (auto const &f : got) {
			if ((f.x - std::floor(f.x) != 0.5f) || (f.y - std::floor(f.y) != 0.5f)) {
				throw Test::error("Rasterizing '" + desc + "', got fragment at (" + std::to_string(f.x) + ", " + std::to_string(f.y) + "), which isn't at a pixel center.");
			}
			if (expected.count(f)) {
				draw(f, '#');
				++matched;
			} else {
				draw(f, '!');
				++extra;
			}
		}
		for (auto const &f : expected) {
			if (!got.count(f)) {
				draw(f, '?');
				++missed;
			}
		}

		if (extra > 0 || missed > 0) {
			//failed!
			std::string info = "Example '" + desc + "' missed " + std::to_string(missed) + " ('?'); had " + std::to_string(extra) + " extra ('!'); and matched " + std::to_string(matched) + " ('#') fragments:";

			//square up the raster:
			size_t width = 0;
			for (auto const &row : raster) {
				width = std::max(width, row.size());
			}
			for (auto &row : raster) {
				row.resize(width, '.');
			}

			for (uint32_t y = static_cast<uint32_t>(raster.size()) - 1; y < static_cast<uint32_t>(raster.size()); --y) {
				info += "\n    " + raster[y];
			}
	
			if (out_of_raster) info += "\n    (" + std::to_string(out_of_raster) + " out-of-range fragments not plotted.)";

			puts(""); //because "test..."
			info("%s", info.c_str());

			throw Test::error("Example '" + desc + "' didn't match expected.");
		}

		//if nothing extra and nothing missed, success!
		assert(matched == expected.size());
	}

	//check that line produces exactly the fragments drawn in a fancy picture:
	void check_line_covers(std::string const &desc, std::initializer_list< Vec2 > const &line_strip, std::initializer_list< std::string > const &raster_) {
		//convert raster to set of points ( with lower-left being (0,0) ):
		std::vector< std::string > raster(raster_);
		std::unordered_set< Vec2 > expected;
		for (uint32_t y = 0; y < raster.size(); ++y) {
			std::string const &row = raster[raster.size()-1-y];
			for (uint32_t x = 0; x < row.size(); ++x) {
				if (row[x] != '.') {
					expected.emplace(x + 0.5f, y + 0.5f);
				}
			}
		}
		//use list-of-points version:
		check_line_covers(desc, line_strip, expected);
	}
}


//--------------------------------------------------

// My test cases
Test test_a1_task2_start_on_diamond_left_point("a1.task2.more.diamond.left", []() {
	extra_test::check_line_covers(
		"line from (3.0f, 1.5f) to (1.25f, 1.5f)",
		{ Vec2(3.0f, 1.5f), Vec2(1.25f, 1.5f) },
		{"....",
		 "..##",
		 "...."}
	);
});

Test test_a1_task2_start_on_diamond_bottom_point("a1.task2.more.diamond.bottom", []() {
	extra_test::check_line_covers(
		"line from (1.5f, 1.0f) to (2.25f, 1.75f)",
		{ Vec2(1.5f, 1.0f), Vec2(2.25f, 1.75f) },
		{"...",
		 ".#.",
		 "..."}
	);
});

Test test_a1_task2_inside_pixel_a_in_b_out("a1.task2.more.inside.pixel.ain.bout", []() {
	extra_test::check_line_covers(
		"line from (1.5f, 1.5f) to (1.2f, 1.2f)",
		{ Vec2(1.5f, 1.5f), Vec2(1.2f, 1.2f) },
		{"...",
		 ".#.",
		 "..."}
	);
});

Test test_a1_task2_test1("a1.task2.more.simple.test1", []() {
 	extra_test::check_line_covers(
  "vertical line",
  { Vec2(1.125f, 1.125f), Vec2(1.125f, 6.875f) },
  {"...",
   ".#.",
   ".#.",
   ".#.",
   ".#.",
   ".#.",
   ".#.",
   "..."}
 );
});

Test test_a1_task2_test2("a1.task2.more.test2", []() {
 	extra_test::check_line_covers(
  "diagonal length 2",
  { Vec2(0.0f, 0.0f), Vec2(2.0f, 2.0f) },
  {".....",
   ".#...",
   "#...."}
 );
});

Test test_a1_task2_test3("a1.task2.more.test3", []() {
 extra_test::check_line_covers(
  "diagonal length 8",
  { Vec2(0.0f, 0.0f), Vec2(8.0f, 8.0f) },
  {".......#.",
   "......#..",
   ".....#...",
   "....#....",
   "...#.....",
   "..#......",
   ".#......",
   "#......."}
 );
}); 

Test test_a1_task2_horizontal_2("a1.task2.more.horizontal.2", []() {
 extra_test::check_line_covers(
  "horizontal line from (1.5, 1.125) to (4.5, 1.125)",
  { Vec2(1.5f, 1.125f), Vec2(4.5f, 1.125f) },
  {".....",
   ".###.",
   "....."}
 );
});

Test test_a1_task2_vertical_2("a1.task2.more.vertical.2", []() {
 extra_test::check_line_covers(
  "vertical line from (1.125, 1.875) to (1.125, 4.125)",
  { Vec2(1.125f, 1.875f), Vec2(1.125f, 4.125f) },
  {"...",
   ".#.",
   ".#.",
   "...",
   "..."}
 );
});

// checks that reaching the top point through the diamond "exits" it
Test test_a1_task2_diamond_top("a1.task2.more.diamond.top", []() {
    extra_test::check_line_covers("line reaches through diamond to its top (1,1)",
                      {Vec2(0.5f, 0.0f), Vec2(1.5f, 2.0f)},
                      {"...", 
						".#.", 
						"#.."});
});

// checks that reaching the top point not through the diamond doesn't "exit" it
Test test_a1_task2_diamond_top_no_exit("a1.task2.more.diamond.top.no.exit", []() {
    extra_test::check_line_covers("line reaches not thorugh diamond to its top (1,1)",
                      {Vec2(0.0f, 1.5f), Vec2(1.5f, 2.0f)},
                      {"...", 
						"#..", 
						"..."});
});

// Checks a horizontal line on the edge of pixels start or end end on the bottom point of diamond
Test test_a1_task2_horizontal_edge_end("a1.task2.more.horizontal.edge.end", []() {
    extra_test::check_line_covers("a horizontal line on the edge of pixels end on the bottom point of diamond",
                      {Vec2(0.8f, 1.0f), Vec2(2.5f, 1.0f)},
                      {"...", 
						".#.", 
						"..."});
});

Test test_a1_task2_horizontal_edge_start("a1.task2.more.horizontal.edge.start", []() {
    extra_test::check_line_covers("a horizontal line on the edge of pixels start on the bottom point of diamond",
                      {Vec2(1.5f, 1.0f), Vec2(2.5f, 1.0f)},
                      {"...", 
						".#.", 
						"..."});
});

// Check 45 degree line

Test test_a1_task2_45("a1.task2.more.45", []() {
    extra_test::check_line_covers("a 45 line from (0.5, 1.0) to (1.5, 2.0)",
                      {Vec2(0.5f, 1.0f), Vec2(1.5f, 2.0f)},
                      {"...", 
						"#.", 
						"..."});
});