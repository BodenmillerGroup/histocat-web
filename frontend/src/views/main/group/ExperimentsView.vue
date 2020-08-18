<template>
  <v-container fluid>
    <v-row class="mt-6 mx-6">
      <v-text-field
        v-model="search"
        append-icon="mdi-magnify"
        label="Search"
        single-line
        hide-details
        clearable
        solo
        class="mt-1"
      />
      <v-select
        v-model="selectedTags"
        :items="tags"
        chips
        deletable-chips
        clearable
        label="Tags"
        multiple
        solo
        class="mt-1 ml-4"
      />
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn v-on="on" dark fab color="primary lighten-1" to="/main/experiments/create" class="ml-4">
            <v-icon>mdi-plus</v-icon>
          </v-btn>
        </template>
        <span>Create experiment</span>
      </v-tooltip>
    </v-row>
    <masonry :cols="{ default: 4, 1000: 3, 700: 2, 400: 1 }" :gutter="{ default: '0px' }">
      <ExperimentCard
        v-for="experiment in experiments"
        :key="experiment.id"
        :experiment="experiment"
        :user="userProfile"
      />
    </masonry>
  </v-container>
</template>

<script lang="ts">
import ExperimentCard from "@/views/main/ExperimentCard.vue";
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";

@Component({
  components: { ExperimentCard },
})
export default class Dashboard extends Vue {
  mainContext = mainModule.context(this.$store);
  experimentContext = experimentModule.context(this.$store);

  search = "";
  selectedTags: string[] = [];

  get tags(): any[] {
    const list = this.experimentContext.getters.tags;
    return list.map((item) => {
      return {
        text: item,
      };
    });
  }

  get userProfile() {
    return this.mainContext.getters.userProfile;
  }

  get experiments() {
    const items = this.search
      ? this.experimentContext.getters.experiments.filter((item) => item.name.includes(this.search))
      : this.experimentContext.getters.experiments;
    if (this.selectedTags.length === 0) {
      return items;
    } else {
      return items.filter((experiment) => {
        if (experiment.tags && this.selectedTags.some((r) => experiment.tags.includes(r))) {
          return experiment;
        }
      });
    }
  }

  async mounted() {
    await Promise.all([this.experimentContext.actions.getExperiments(), this.experimentContext.actions.getTags()]);
  }
}
</script>
