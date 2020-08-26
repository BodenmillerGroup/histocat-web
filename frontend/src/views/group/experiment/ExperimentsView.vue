<template>
  <v-col>
    <v-toolbar dense class="toolbar">
      <v-toolbar-title>Experiments</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn
          v-if="isGroupAdmin"
          text
          color="primary"
          :to="{ name: 'group-experiments-create', params: { groupId: activeGroupId } }"
        >
          Add Experiment
        </v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-expansion-panels>
      <v-expansion-panel>
        <v-expansion-panel-header>Filter</v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-select
            v-model="tagsFilter"
            :items="tags"
            item-text="text"
            item-value="value"
            chips
            clearable
            label="Tags"
            multiple
            prepend-icon="mdi-filter-outline"
            solo
            dense
          >
            <template v-slot:selection="{ attrs, item, select, selected }">
              <v-chip v-bind="attrs" :input-value="selected" close @click="select" @click:close="removeTagFilter(item)">
                {{ item }}
              </v-chip>
            </template>
          </v-select>
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
    <v-card>
      <v-card-title>
        <v-spacer />
        <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="items"
        :loading="!items"
        :search="search"
        :custom-filter="filter"
        :items-per-page="15"
        :footer-props="{
          itemsPerPageOptions: [10, 15, 20, -1],
          showFirstLastPage: true,
          showCurrentPage: true,
        }"
        multi-sort
      >
        <template v-slot:item.name="{ item }">
          <router-link class="link" :to="{ name: 'group-experiment', params: { experimentId: item.id } }">
            {{ item.name }}
          </router-link>
        </template>
        <template v-slot:item.tags="{ item }">
          <v-chip v-for="tag in item.tags" :key="tag" x-small pill disabled class="mr-1">
            {{ tag }}
          </v-chip>
        </template>
        <template v-slot:item.created_at="{ item }">
          {{ item.created_at | stringToUTCString }}
        </template>
        <template v-slot:item.action="{ item }">
          <v-menu bottom left>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on">
                <v-icon>mdi-dots-vertical</v-icon>
              </v-btn>
            </template>
            <v-list dense>
              <v-list-item :to="{ name: 'group-experiments-edit', params: { experimentId: item.id } }">
                <v-list-item-icon>
                  <v-icon color="grey">mdi-pencil-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Edit</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item v-if="isGroupAdmin" @click="deleteExperiment(item.id)">
                <v-list-item-icon>
                  <v-icon color="red accent-1">mdi-delete-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Delete</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list>
          </v-menu>
        </template>
      </v-data-table>
    </v-card>
  </v-col>
</template>

<script lang="ts">
import { experimentModule } from "@/modules/experiment";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

@Component
export default class Dashboard extends Vue {
  mainContext = mainModule.context(this.$store);
  groupContext = groupModule.context(this.$store);
  experimentContext = experimentModule.context(this.$store);

  search = "";
  tagsFilter: string[] = [];

  readonly headers = [
    {
      text: "Name",
      value: "name",
    },
    {
      text: "Description",
      value: "description",
    },
    {
      text: "Tags",
      value: "tags",
      filterable: false,
      sortable: false,
    },
    {
      text: "Created",
      value: "created_at",
      filterable: false,
      width: "240",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
      filterable: false,
      width: "70",
    },
  ];

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

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

  get items() {
    let items = this.experimentContext.getters.experiments;
    if (this.tagsFilter.length > 0) {
      items = items.filter((experiment) => {
        if (experiment.tags && this.tagsFilter.some((r) => experiment.tags.includes(r))) {
          return experiment;
        }
      });
    }
    return items;
  }

  filter(value, search, item) {
    if (!search) {
      return true;
    }
    const normalizedSearchTerm = search.toLowerCase().trim();
    return item.name.toLowerCase().indexOf(normalizedSearchTerm) !== -1 || item.description
      ? item.description.toLowerCase().indexOf(normalizedSearchTerm) !== -1
      : false;
  }

  removeTagFilter(item) {
    this.tagsFilter.splice(this.tagsFilter.indexOf(item), 1);
    this.tagsFilter = [...this.tagsFilter];
  }

  async deleteExperiment(id: number) {
    if (self.confirm("Are you sure you want to delete the experiment?")) {
      await this.experimentContext.actions.deleteExperiment(id);
    }
  }

  async mounted() {
    await Promise.all([
      this.experimentContext.actions.getGroupExperiments(this.activeGroupId!),
      this.experimentContext.actions.getExperimentTags(this.activeGroupId!),
    ]);
  }
}
</script>

<style scoped>
.toolbar {
  margin-bottom: 10px;
}
</style>
